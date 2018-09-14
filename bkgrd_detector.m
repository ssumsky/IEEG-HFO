function bkgrd_detector(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: The function processes a *.mat data file (name: "filename") provided by the
%              user to extract the onset and offset time of putative HFO events in the
%              common reference signal for sake of artifact removal. It is assumed that
%              the file contains a multichannel time series (double) with time labels. The
%              HFOs are determined by using the Staba algorithm (Staba et al., J. Neuro-
%              physiol., 2002). The onset and offset time of each HFO are returned in 2-D
%              MATLAB arrays (one array per sequence) and stored in a *.mat file in the
%              current folder.
%
%-----------------------------------------------------------------------------------------
%
% Input:    filename - String of characters reporting the name of the *.mat file.
%
%-----------------------------------------------------------------------------------------
%
% Output:   no output variable returned.
%
%-----------------------------------------------------------------------------------------
%
% List of notable invoked MATLAB functions:
%
%               a) dsp.BiquadFilter()     [ from the Signal Processing Toolbox ]
%               b) findpeaks()            [ from the Signal Processing Toolbox ]
%               c) iirnotch()             [ from the DSP System Toolbox ]
%
%-----------------------------------------------------------------------------------------
%
% Author: S. Santaniello
%
% Ver.: 3.0 - Date: 08/03/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%----------------------------------------------------------------------------------------%
% 1. Check compulsory input variables and load required information
%----------------------------------------------------------------------------------------%
%
% check the first input variable that is provided
if (nargin<1)
    error('FILENAME OF THE TIME SERIES MUST BE PROVIDED'); 
else
    if (~ischar(filename))
        error('FILENAME MUST BE A STRING'); 
    else
        if (~contains(filename,'.mat'))
            error('FILENAME MUST BE VALID');
        end
    end
end

% load the data. Note:
%   - "data"   -> it stores the iEEG time series (one time series per column);
%   - "time"   -> it stores the time marks (in seconds) of the iEEG samples;
%   - "labels" -> it stores the name of each iEEG channels (one label per channel);
load(filename,'-regexp','data|time|labels');
%----------------------------------------------------------------------------------------%



%----------------------------------------------------------------------------------------%
% 2. Initialization and definition of environment variables
%----------------------------------------------------------------------------------------%
%
% sampling rate (Hz)
FS = round(1/(time(2)-time(1)));

% size of each window (in number of samples)
WIN = round(600*FS);

% size of the shift between consecutive windows (in number of samples)
SHFT = WIN;

% size of the pre-/post-window intervals( in number of samples)
DELT = round(0.25*FS);

% Notch filters to remove the power-line component and its secondary harmonics
BW = 3;                        % bandwidth (Hz)
Apass = 1;                     % bandwidth attenuation (in dB)
Fnotch = 60:60:(FS/2);         % Notch frequencies (in Hz)
for i=1:length(Fnotch)
    eval(sprintf('[num%d,den%d] = iirnotch(%f,%f,%f);',i,i,Fnotch(i)/(FS/2),BW/(FS/2),Apass));
end

% Band-pass filter with cutoff at 100 Hz and 500 Hz
% % Fstop1 = 20;    % First Stopband Frequency
% % Fpass1 = 100;   % First Passband Frequency
% % Fpass2 = 500;   % Second Passband Frequency
% % Fstop2 = 1000;  % Second Stopband Frequency
% % Astop1 = 60;    % First Stopband Attenuation (dB)
% % Apass  = 1;     % Passband Ripple (dB)
% % Astop2 = 80;    % Second Stopband Attenuation (dB)
% % h      = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass, Astop2, FS);
% % Hd     = design(h, 'butter', 'MatchExactly', 'stopband', 'SystemObject', true);
Nord  = 6;      % Order
F3dB1 = 100;    % First
F3dB2 = 500;    % Second
h     = fdesign.bandpass('n,f3db1,f3db2', Nord, F3dB1, F3dB2, FS);
Hd    = design(h, 'butter', 'SystemObject', true);

% parameters to be used for HFO detection:
%   a) nrms   = number of samples used to compute the RMS;
%   b) mindur = minimum duration (in no. of samples) of a HFO;
%   c) mindis = minimum separation (in no. of samples) between two consecutive HFOs;
%   d) maxdur = maximum duration (in no. of samples) of a HFO;
%   e) npks   = minimum number of supra-threshold peaks in the rectified HFO.
nrms   = round(0.003*FS); mindur = round(0.006*FS);
mindis = round(0.010*FS); npks   = 6;
%-----------------------------------------------------------------------------------------



%-----------------------------------------------------------------------------------------
% 3. Extraction of the HFOs
%-----------------------------------------------------------------------------------------
%
% onset and end time (in sec) of the series. It is required to account for gaps in the
% series 
st_pnt = time(1); end_pnt = time(end);

% enlarge the data set with constant values at the onset and end of data series to remove
% artifactual oscillations
data = [repmat(data(1,:),DELT,1); data; repmat(data(end,:),DELT,1)];
time = [(-DELT:-1)'.*(1/FS)+time(1); time; time(end)+(1/FS).*(1:DELT)'];

% apply a moving median filter to remove artifactual spikes
data = movmedian(data,5,1);

% initialize the array for the HFO events
HFO_bkgrd = zeros(0,2);

% loop across the consecutive windows...
tic
while (st_pnt<end_pnt)
    
    % set the pointers to the data samples collected within the interval of interest
    pntval  = (time>=st_pnt-DELT*(1/FS) & time<st_pnt+WIN*(1/FS));
    pntval0 = (time>=st_pnt & time<st_pnt+WIN*(1/FS));
    mask    = pntval0(pntval);
    
    % extract the data of interest
    vls = data(pntval,:);
    tvl = time(pntval0);
    
    % check that in the interval at least some channels are connected
    flg = (sum(sum(vls==0,1)==size(vls,1))<size(vls,2));
    
    % if there is enough points in the interval...
    if (sum(pntval0)>0.4*WIN && flg)
            
        % identify the channels that were disconnected during the recording
        pos = (sum(vls==0,1)==size(vls,1));
        
        % compute the common average reference
        ref = mean(vls(:,~pos),2);
        
        % pre-filter the common reference signal and remove the pre-/post-interval
        % data points
        for i=1:length(Fnotch)
            ref = filtfilt(eval(sprintf('num%d',i)),eval(sprintf('den%d',i)),ref);
        end
        ref = filtfilt(Hd.SOSMatrix,Hd.ScaleValues,ref);
        vls = ref(mask);
        clear pos ref
        
        % compute the RMS amplitude on consecutive sub-windows
        rmsval = vls.^2;
        for w=2:nrms
            rmsval = rmsval + [zeros(w-1,1); vls(1:end-w+1,1).^2];
        end
        rmsval = sqrt(rmsval./nrms);
        
        % compute the rectified signals
        absval = abs(vls);
        
        % compute mean and SD for the RMS and rectified ECoG time series
        avg_rms = mean(rmsval); std_rms = std(rmsval,[]);
        avg_abs = mean(absval); std_abs = std(absval,[]);
        
        
        % A. detect candidate HFO artifacts. Criteria include:
        %---------------------------------------------------------------------------------
        %   a) RMS > MEAN + 5*SD;
        %   b) condition a) valid for "mindur" consecutive samples;
        %   c) HFOs must be at least "mindis" samples apart.
        %
        % Candidate events are stored in a 2-D array along with the info:
        %   - onset point;
        %   - offset (i.e., end) point.
        candidate_HFO = zeros(0,2);
        
        % check condition a)
        tmp1 = (rmsval - avg_rms)./std_rms;
        tmp1(isnan(tmp1)) = 0;
        pnt = diff([0; (tmp1 > 5)]);
        
        % pointers to the onset and offset of each candidate HFO
        onset = find(pnt==1);
        offst = find(pnt==-1);
        if (length(offst)<length(onset))
            eval(sprintf('offst(end+1,1) = %d;',size(tmp1,1)));
        end
        
        % check condition b)
        if (~isempty(onset))
            pos = ((offst - onset) > mindur);
            onset = onset(pos);
            offst = offst(pos);
            clear pos
        end
        
        % check condition c)
        if(~isempty(onset))
            pos = ((onset(2:end) - offst(1:end-1))<mindis);
            onset = onset([true; ~pos]);
            offst = offst([~pos; true]);
            clear pos
        end
        
        % store the candidate HFO events
        candidate_HFO(end+1:end+length(onset),:) = [onset(:) offst(:)];
        clear onset offst tmp1 pnt rmsval avg_rms std_rms
        %---------------------------------------------------------------------------------
        
        
        % B. detect likely HFO artifacts. Criteria include:
        %---------------------------------------------------------------------------------
        %   d) rectified HFO must include at least 3 peaks;
        %   e) peaks of rectified HFO must be > MEAN + 3*SD.
        %
        % HFO events are stored in a 2-D array along with the info:
        %   - onset time  (in microseconds);
        %   - offset time (in microseconds).
        final_HFO = zeros(0,2);
        
        % check conditions d)
        [~,locs] = findpeaks(absval,'MinPeakHeight',avg_abs+3*std_abs);
        
        % check condition e)
        if (~isempty(locs))
            
            % for each candidate HFO artifact in the list...
            for r=1:size(candidate_HFO,1)
                lcsv = sum(locs>=candidate_HFO(r,1) & locs<=candidate_HFO(r,2));
                
                % save the identified HFO event
                if (lcsv>=npks)
                    final_HFO(end+1,:) = [tvl(candidate_HFO(r,1)) tvl(candidate_HFO(r,2))];
                end
                clear lcsv
            end
        end
        clear locs absval avg_abs std_abs
        %---------------------------------------------------------------------------------
        
        % save the results
        HFO_bkgrd(end+1:end+size(final_HFO,1),:) = final_HFO;
        
        % update the pointer to the onset of next window of data
        st_pnt = st_pnt+SHFT*(1/FS);
    else
        
        % update the pointer to the next available data sample
        st_pnt = time(find(time>=st_pnt+WIN*(1/FS),1,'first'));
    end
    clear pntval pntval0 vls tvl
end

% save the results in repository *.mat files
if ~exist(sprintf('staba_times_%s',filename),'file')
    save(sprintf('staba_times_%s',filename),'-regexp','^HFO_');
else
    save(sprintf('staba_times_%s',filename),'-regexp','^HFO_','-append');
end

elapsed = toc;
fprintf('%s: elapsed time %f\n',filename,elapsed);
%-----------------------------------------------------------------------------------------
