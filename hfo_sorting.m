function hfo_sorting(filename)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: The function processes HFO events stored in the *.mat file "filename" and
%              sorts these events into different arrays based on channels and existing
%              holes in the original iEEG multichannel time series. Results are stored in
%              the same file.
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
% List of functions that must have been invoked beforehand:
%
%               a) bkgrd_detector()     [ from the current folder ]
%               b) staba_detector()     [ from the current folder ]
%               c) popDet_detector()    [ from the current folder ]
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
    error('FILENAME OF THE HFO TIME SERIES MUST BE PROVIDED'); 
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
%   - "FS"        -> iEEG sampling rate (Hz);
%   - "ncol"      -> number of iEEG channels;
%   - "HFO"       -> Nx3 matrix that stores the candidate HFOs [onset, offset, channel #];
%   - "aHFO"      -> Nx3 matrix that stores HFOs >800 Hz [onset, offset, channel #];
%   - "HFO_bkgrd" -> Mx2 matrix that stores common reference HFOs [onset, offset];
load(filename,'-regexp','HFO|aHFO|HFO_bkgrd|ncol');
%----------------------------------------------------------------------------------------%



%----------------------------------------------------------------------------------------%
% 2. Initialization and definition of environment variables
%----------------------------------------------------------------------------------------%
%
% initialize the struct array for the HFO events and counts. The structure includes one
% entry per channel and each entry has one or more fields labelled "series" (one for each
% sequence of continuous data). Each "series" includes three fields:
%
%   - marks: time marks (in seconds) of onset and termination of each detected HFO;
%   - count: Nx1 vector including the number of HFOs detected in N consecutive intervals,
%            each interval being of size WIN;
%   - time : Nx1 vector of time marks (in seconds) corresponding to the end of each
%            interval considered in "count".
HFOevents(ncol) = struct('series',[]);
for k=1:ncol
    HFOevents(k).series(1).count = zeros(0,1);
    HFOevents(k).series(1).marks = zeros(0,2);
end
%-----------------------------------------------------------------------------------------



%-----------------------------------------------------------------------------------------
% 3. Extraction of the HFOs
%-----------------------------------------------------------------------------------------
%    
% Step 1) Remove HFOs that are likely artifacts (i.e., within 100ms from any HFO detected
% in the common reference signal). Note that time values are reported in seconds
termflg = false; k = 0; count = 0;
while (k<size(HFO_bkgrd,1) && ~termflg)
    if (~isempty(HFO(:)))
        k = k+1;
        pos = (HFO(:,2)<(HFO_bkgrd(k,1)-0.1) | HFO(:,1)>(HFO_bkgrd(k,1)+0.1));
        HFO = HFO(pos,:);
        fprintf('Bkgrd HFO %d... %d HFOs pruned\n',k,sum(~pos));
        count = count + sum(~pos);
    else
        warning('No HFO detected in the block');
        termflg = true;
    end
end
if (isempty(HFO(:)))
    fprintf('%s: No HFO in the file after pruning...\t',filename);
    termflg = true;
end

% Step 2) Remove HFOs that are likely artifacts (i.e., within 100ms from any HFO detected
% beyond 800 Hz). Note that time values are reported in seconds
termflg2 = false;
if (~termflg)
    k = 0;
    while (k<size(aHFO,1) && ~termflg2)
        if (~isempty(HFO(:)))
            k = k+1;
            pos = ((HFO(:,2)<(aHFO(k,1)-0.1) | HFO(:,1)>(aHFO(k,1)+0.1)) & HFO(:,3)==aHFO(k,3)) | HFO(:,3)~=aHFO(k,3);
            HFO = HFO(pos,:);
            fprintf('aHFO %d... %d HFOs pruned in channel %d\n',k,sum(~pos),aHFO(k,3));
            count = count + sum(~pos);
        else
            warning('No HFO detected in the block');
            termflg2 = true;
        end
    end
    if (isempty(HFO(:)))
        fprintf('%s: No HFO in the file after pruning...\t',filename);
        termflg2 = true;
    end
end

% Step 3) Extract the HFO-related info and count the number of HFOs per channel in each
% interval. If an interval was NOT processed because of lack of data, break the count
% series
if (~termflg && ~termflg2)

    % count the number of HFOs in each channel in the current interval
    for k=1:ncol
        
        % pointer to the events occurring in the current window
        tmpHFO = (HFO(:,3)==k);
        
        % save the detected events
        HFOevents(k).series(1).count(end+1,1) = sum(tmpHFO);
        if (sum(tmpHFO)>0)
            HFOevents(k).series(1).marks(end+1:end+sum(tmpHFO),:) = HFO(tmpHFO,1:2);
        end
        clear tmpHFO
    end
end
fprintf('Total: %d HFOs pruned across channels\n',count);

% save the results in repository *.mat files
save(filename,'-regexp','^HFOevents','-append');
%-----------------------------------------------------------------------------------------
