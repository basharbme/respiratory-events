function [ apnea_events_epochs,extrema,signal_resampled,apnea_events_resampled,rsf] = sn_findRespiratoryEvents(varargin)
%reads varargins of a function and gives back the parsed parameters Compumedics dpsg files and converts in matlab struct
%
% cli:
%   cwlVersion: v1.0-extended
%   class: matlabfunction
%   baseCommand: [events,extrema] = sn_findRespiratoryEvents(varargin)
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
%     lowerlimitmovingaverageamplitude:
%       type: float?
%       inputBinding:
%         prefix: llma
%       doc: "The lower limit for accepted amplitudes (Peak2Peak) with respect to moving
%         average of amplitudes, as fraction of the moving average. Default 0.5"
%     apneaminimumduration:
%       type: float?
%       inputBinding:
%         prefix: amd
%       doc: "The lower limit for apnea duration in seconds. 
%               We measure durations fthe exact timepoints as in AASM, so the value
%               should be higher than 10 secs. Default 15 secs"
%     debug:
%       type: int?
%       inputBinding:
%         prefix: debug
%       doc: "if set to 1 debug information is provided. Default 0"
%
%   outputs:
%     apnea_events_epochs:
%       type: matlabArray
%       doc: "An array with start and stop epochs of a detected apnea."
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
% Notes
% The criteria along AASM are the following:
%  apnea: <=90% peak amplitude for >= 10 seconds relative to pre-event
%  baseline
%  hypopnea: <=30% peak amplitude for >= 10 seconds relative to pre-event
%  baseline AND >= 3% drop in Sa02 or arousal is assosiated 
%  we take only the airflow, so hypopneas can only be candidates, that must
%  be checked against Sa02

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 0. Parse Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% required input
myinput.data = NaN;
% sample frequency
myinput.sf = 1;
% resample frequency
myinput.rsf = 4;
% lower limit accepted amplitude as breathing
myinput.llma = 0.1;
% minimum period
myinput.mp = 2;
% apnea minimum duration
myinput.amd = 15;
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
%% 2. Get amplitudes and durations of breathing cycles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find extreme values with constraints on the minimum distance in time between two maxima. 
[breathingrate,extrema,signal_resampled]  ...
    = sn_getRespiratoryRate('data',myinput.data,'sf',myinput.sf,...
                            'llma',myinput.llma,'rsf',myinput.rsf,...
                            'mp',myinput.mp);
                        
% get the duration of the found breathing cycles, 
% as the duration is given in samples of the signal_resampled, it needs to
% be converted in seconds
duration_bc = diff(extrema(:,1))/myinput.rsf;

%find durations higher than the minimum duration set
apneas = duration_bc >= myinput.amd;
if (myinput.debug), whos apneas, end

%get the start and end in seconds
apnea_events_resampled = [extrema([apneas;logical(0)],1),extrema([logical(0);apneas],1)];

%get the start and end in terms of epochs
apnea_events_epochs = fix(apnea_events_resampled/(myinput.rsf*30))+1;

%output rfs
rsf = myinput.rsf;
    
    
