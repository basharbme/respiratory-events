function [breathingrate,extrema,signal_resampled]  = sn_getRespiratoryRate(varargin)
%calculates the respiratory rate - its the successor of sn_getBreathingRate, that is kept for downwards compatibility 
%
% cli:
%   cwlVersion: v1.0-extended
%   class: matlabfunction
%   baseCommand: [breathingrate,extrema,signal_resampled] = sn_getRespiratoryRate(varargin)
%
%   inputs:
%     data:
%       type: matlabArray
%       inputBinding:
%         prefix: data
%       doc: "A matlab array containing the airflow."
%     samplingfrequency:
%       type: float?
%       inputBinding:
%         prefix: sf
%       doc: "The sampling frequency in Hz. Default: 1 Hz"
%     minimumperiod:
%       type: float?
%       inputBinding:
%         prefix: mp
%       doc: "The minimum period of signal. Default: 2 secs"
%     resamplingfrequency:
%       type: float?
%       inputBinding:
%         prefix: rsf
%       doc: "The resampling frequency. Default: 4 Hz"
%     numbersamplesmovingaverage:
%       type: float?
%       inputBinding:
%         prefix: nsma
%       doc: "The number of samples in moving average. Default 11 samples"
%     lowerlimitmovingaverageduration:
%       type: float?
%       inputBinding:
%         prefix: llmd
%       doc: "The lower limit for accepted durations (Peak2Peak) with respect to moving
%         average of duration, as fraction of the moving average. Default 0.5"
%     lowerlimitmovingaverageamplitude:
%       type: float?
%       inputBinding:
%         prefix: llma
%       doc: "The lower limit for accepted amplitudes (Peak2Peak) with respect to moving
%         average of amplitudes, as fraction of the moving average. Default 0.5"
%     breathingratesamplingfrequency:
%       type: float?
%       inputBinding:
%         prefix: brsf
%       doc: "The sampling frequency of breathing rate. Default 1 Hz"
%     debug:
%       type: int?
%       inputBinding:
%         prefix: debug
%       doc: "if set to 1 debug information is provided. Default 0"
%
%   outputs:
%     breathingrate:
%       type: matlabArray
%       doc: "vector containing the breathing rate"
%     extrema:
%       type: matlabArray
%       doc: "matrix containing the found extreme values
%               cols: 
%               1: location maximum (in resampled signal)         
%               2: value maximum         
%               3: location minimum (in resampled signal)         
%               4: value minimum"         
%
%
%   s:author:
%     - class: s:Person
%       s:identifier:  https://orcid.org/0000-0002-7238-5339
%       s:email: mailto:dagmar.krefting@htw-berlin.de
%       s:name: Dagmar Krefting
% 
%   s:dateCreated: "2018-10-03"
%   s:license: https://spdx.org/licenses/Apache-2.0 
% 
%   s:keywords: edam:topic_3063, edam:topic_2082
%     doc: 3063: medical informatics, 2082: matrix
%   s:programmingLanguage: matlab
% 
%   $namespaces:
%     s: https://schema.org/
%     edam: http://edamontology.org/
% 
%   $schemas:
%     - https://schema.org/docs/schema_org_rdfa.html
%     - http://edamontology.org/EDAM_1.18.owl
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Parse Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% required input
myinput.data = NaN;
% sample frequency
myinput.sf = 1;
% resample frequency
myinput.rsf = 4;
% minimu period
myinput.mp = 2;
% number of samples in moving average
myinput.nsma = 11;
% lower limit of amplitude in moving average
myinput.llma = 0.5;
% lower limit of duration in moving average
myinput.llmd = 0.5;
% breathing rate sampling frequency
myinput.brsf = 1;
% debug level off
myinput.debug = 0;

try
    myinput = mt_parameterparser('myinputstruct',myinput,'varargins',varargin);
catch ME
    disp(ME)
    return
end

if (myinput.debug)
    myinput
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Resample signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% resample the signal to rsf
signal_resampled = resample(myinput.data,myinput.rsf,myinput.sf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Find extreme values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find extreme values with constraints on the minimum distance in time between two maxima. 
% It is always the minimum preceding the maximum in one row. 

extrema = sn_getExtrema(signal_resampled,'mp',myinput.mp,'sf',myinput.rsf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Exclude false extreme values with small amplitudes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% how we define "false" breathing cycles?

%% Check for amplitudes (peak2peak)

%get distance between maximum and minimum (amplitude, peak2peak)
extrema_dist_min_max = extrema(:,2)-extrema(:,4);
%now do the same for the amplitude between maximum and following minimum
%this array has one element less, as it compares with the following row,
%and the last row does not have a following one
extrema_dist_max_min = extrema(1:end-1,2)-extrema(2:end,4);

%get moving average of amplitudes (peak2peak)
% we use only the first amplitudes
extrema_dist_mavg = mt_movingAverage('data',extrema_dist_min_max,'ns',myinput.nsma);

%get extrema with amplitudes smaller then lower limit 
extrema_dist_bin = extrema_dist_min_max < myinput.llma*extrema_dist_mavg;

%% Check for breathing cycle durations
% This is relevant for detection of respiratory events rather than for the
% breathing rate, as apneas would be seen as no breathing, leading to a
% very low breathing rate, and that's okay. 

%we check this only for the maxima in a first attempt
duration_dist_max = diff(extrema(:,1));
%get moving average
duration_dist_max_mavg = mt_movingAverage('data',duration_dist_max,'ns',myinput.nsma);
%get extrema with durations smaller then lower limit 
duration_dist_bin = duration_dist_max < myinput.llmd*duration_dist_max_mavg;

%% Check for absolute values with respect to preceding maximum and following minimum

%get extrema with larger maximum than preceding extremum
extrema_max_ltp_bin = extrema(:,2) > [0; extrema(1:end-1,2)];

%get extrema with smaller minimum than following extremum
extrema_min_ltf_bin = extrema(:,4) < [extrema(2:end,4); 0];

%% before deleting false extrema, correct for larger maxima and smaller minima

%false positives with larger maxima (but small amplitude!)
fplm = (extrema_dist_bin & extrema_max_ltp_bin);

%false positives with smaller minima (but small amplitude!)
fpsm = (extrema_dist_bin & extrema_min_ltf_bin);

%%  DEBUGGING
%find those, where both conditions are given, they must be excluded from
%shifting, otherwise extrema order is not preserved
%probably there is a better solution, but I don't see them in the moment.
%Given the rare situation, the resulting error might be okay

%% find extrema, where absolute values are larger than neighbors

fpsm_fplm = fpsm & fplm;
%invert to set incidents to false
fpsm_fplm = ~fpsm_fplm;
%remove these from fpsm and fplm
fpsm = fpsm & fpsm_fplm;
fplm = fplm & fpsm_fplm;

%% debugging finished

%store the maximal value to precessor (keep the extreme value)
extrema([fplm(2:end); logical(0)],1:2) = extrema([logical(0); fplm(2:end)],1:2);

%store the minimal value to follower (keep the extreme value)
extrema([logical(0); fpsm(1:end-1)],3:4) = extrema([fpsm(1:end-1); logical(0)],3:4);

% delete extreme values below fraction of moving average
extrema(extrema_dist_bin,:) = [];

%delete also the moving average-values
extrema_dist_mavg(extrema_dist_bin,:) = [];

%% and now the same story for the amplitudes between preceding maxima and following minima 

%get distance between maximum and following minimum
extrema_dist = extrema(1:end-1,2)-extrema(2:end,4);

%use the same moving average as before!!!
%get extrema with distances smaller then lower limit 
edb = extrema_dist < myinput.llma*extrema_dist_mavg(1:end-1);

%write minimum values from these false positives to follower
extrema([logical(0); edb(1:end-1)],3:4) = extrema([edb(1:end-1); logical(0)],3:4);

%delete false positives 
extrema(edb,:) = [];

%% Get breathing rate

%calculate period: diff between minima and maxima
bp = diff(extrema(:,[1,3]))/myinput.rsf;

%moments to which the calculated period should be assigned: middle point
%between the two time points the period is calculated from
bpi = extrema(1:end-1,[1,3])+bp/2;

%flip cols for having the preceding minimum earlier than following maximum
bp = reshape(flipdim(bp,2)',1,2*length(bp));
bpi = reshape(flipdim(bpi,2)',1,2*length(bpi));
whos bpi
%interpolate breathingperiods between first and last extrem value;
br = interp1(bpi,(1./bp),(floor(bpi(1)):ceil(bpi(end))));

%% needs some correction
%pad at both ends to fit resampled data
padstart = floor(bpi(1))-1;
padend = length(signal_resampled)-ceil(bpi(end));
%if floor is 0, cut first point, not to get confused with indices and
%timepoints
if (padstart == -1)
breathingrate = [br(2:end),repmat(br(end),1,padend)];
else    
breathingrate = [repmat(br(1),1,padstart),br,repmat(br(end),1,padend)];
end

%resample breatingrate to brsf
breathingrate = resample(breathingrate,myinput.brsf,myinput.rsf);

end






