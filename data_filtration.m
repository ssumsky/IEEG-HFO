function data_filtration(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: The function processes a *.mat data file (name: "filename") provided by the
%              user by applying notch- and band-pass (100-500 Hz) filters to each channel.
%              It is assumed that the file contains a multichannel time series (double)
%              with time labels. The filtered time series are stored in a *.mat file in
%              the current folder.
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
%               b) iirnotch()             [ from the DSP System Toolbox ]
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
% number of recording channels
ncol = size(data,2);

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

% Band-pass filter with cutoff at 100 Hz and 500 Hz (8th order Butterworth filter)
Hd = dsp.BiquadFilter('Structure', 'Direct form II', ...
    'SOSMatrix', [1 0 -1 1 -1.77101600844797 0.843504967673282; 1 0 -1 1 ...
    0.170951124201346 0.533946648171601; 1 0 -1 1 -1.46331985630811 ...
    0.55370082215595; 1 0 -1 1 -0.0956099995668841 0.0803124533188455], ...
    'ScaleValues', [0.558588907862343; 0.558588907862343; 0.474021677384833; 0.474021677384833; 1]);
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
data = medfilt1(data,5,[],1);

data = [repmat(data(1,:),DELT,1); data; repmat(data(end,:),DELT,1)];
time = [(-DELT:-1)'.*(1/FS)+time(1); time; time(end)+(1/FS).*(1:DELT)'];

% initialize the array for the filtered data
data_filt = zeros(size(data));

% loop across the consecutive windows...
tic
while (st_pnt<end_pnt)
    
    % set the pointers to the data samples collected within the interval of interest
    pntval  = (time>=st_pnt-DELT*(1/FS) & time<st_pnt+(WIN+DELT)*(1/FS));
    pntval0 = (time>=st_pnt & time<st_pnt+WIN*(1/FS));
    mask    = pntval0(pntval);
        
    % extract the data of interest
    vls = data(pntval,:);
        
    % check that in the interval at least some channels are connected
    flg = (sum(sum(vls==0,1)==size(vls,1))<size(vls,2));
    
    % if there is enough points in the interval...
    if (sum(pntval0)>0.4*WIN && flg)

        % identify the channels that were disconnected during the recording
        pos = (sum(vls==0,1)==size(vls,1));
        
        % apply the common average reference
        ref = mean(vls(:,~pos),2);
        vls(:,~pos) = vls(:,~pos) - repmat(ref,1,ncol-sum(pos));
        
        % pre-filter the time series and remove the pre-/post-interval data points
        for i=1:length(Fnotch)
            vls(:,~pos) = filtfilt(eval(sprintf('num%d',i)),eval(sprintf('den%d',i)),vls(:,~pos));
        end
        vls = step(Hd,vls); vls(:,pos) = 0;
        vls = vls(mask,:);
        data_filt(pntval0,:) = vls;
        clear pos ref
        
        % update the pointer to the onset of next window of data
        st_pnt = st_pnt+SHFT*(1/FS);
    else
        
        % update the pointer to the next available data sample
        st_pnt = time(find(time>=st_pnt+WIN*(1/FS),1,'first'));
    end
    clear pntval pntval0 vls tvl
end

% resize the filtered data series
data_filt = data_filt(DELT+1:end-DELT,:);

% save the results in repository *.mat files
if ~exist(sprintf('filt_%s',filename),'file')
    save(sprintf('filt_%s',filename),'-regexp','^data_filt');
else
    save(sprintf('filt_%s',filename),'-regexp','^data_filt','-append');
end

elapsed = toc;
fprintf('%s: elapsed time %f\n',filename,elapsed);
%-----------------------------------------------------------------------------------------
