function [PAFrames, flw, USFrames, pathname] = LoadFiles(varargin)
%function settings
nr_load_fr=2; %how many consectives frames in file
nr_skip_fr=18; %how many frames to skip reading
ld_existing = 0;

if ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'Path'
                pathname = varargin{input_index + 1};
            case 'LoadExisting'
                ld_existing = varargin{input_index + 1};
            otherwise
                error('Unknown optional input');
        end
    end
end              
%select folder
if ~exist('pathname','var')
    pathname = uigetdir('/Users/Thore/Documents/LabData/Transducer Measurements/160707LinArray/');
%     pathname = uigetdir('C:\Users\LABPC_TB\Documents\TransducerMeasurements\201605');
end
%% Initialise
%create list of folder contents
cases=struct2cell(dir(pathname));
cases=cases(1,3:end);

%load folder settings
[flow_settings, flow_rate, RF_settings] = LoadSettings(pathname);

%Find PA and US signals
PAFiles={};USFiles={};
for i=1:length(cases)
    if strfind(cases{i},'flow_PA_');
        PAFiles=[PAFiles; cases{i}];
    elseif strfind(cases{i},'flow_US_');
        USFiles=[USFiles; cases{i}];
    end
end

%% Load data
% load PA file settings
PAFrames = LoadFrames(pathname, 'PAFrames', PAFiles,...
    nr_load_fr, nr_skip_fr,flow_settings, flow_rate, RF_settings, ld_existing);
if ~isa(PAFrames,'Frames')
    disp('changing type')
    PAFrames = typecast(PAFrames,'Frames');
end

% load US file settings
USFrames = LoadFrames(pathname, 'USFrames', USFiles,...
    20,0,flow_settings, flow_rate, RF_settings, ld_existing);
if ~isa(USFrames,'Frames')
    disp('changing type')
    USFrames = typecast(USFrames,'Frames');
end


PAFrames.KWaveInit();
flw.theta=PAFrames.flw.theta;
flw.tube_diameter=PAFrames.flw.tube_diameter;
flw.PRF=PAFrames.flw.PRF;
flw.rate=PAFrames.flw_r;
USFrames.KWaveInit();
temp = RFObj(PAFrames.acq.fs,PAFrames.RF.speed_of_sound,...
    PAFrames.RF.pitch*1E-3,PAFrames.pathname,'RFobj',...
    size(PAFrames.rfm,2),size(PAFrames.rfm,1));
temp.rfm=PAFrames.rfm;
PAFrames=temp;
end

function [frames] = LoadFrames(pathname,objname,files,...
    nr_load_fr, nr_skip_fr,flow_settings, flow_rate, RF_settings, ld_existing)
currentpath=pwd;
cd(pathname)
if exist([objname,'.mat'], 'file') && ~ld_existing
    title=['Load ',objname];
    msg=['Load existing ',objname,'.mat?'];
    button=questdlg(msg,title,'Yes','No','Yes');
    switch button
        case 'Yes'
            frames = load([pathname,'/',objname,'.mat']);
            cd(currentpath)
            frames = frames.obj;            
            return
        case 'No'
            % load file settings
    end
elseif exist([objname,'.mat'], 'file') && ld_existing
    frames = load([pathname,'/',objname,'.mat']);
    cd(currentpath)
    frames = frames.obj;
    return
end
cd(currentpath)
 % load file settings
    filename=[pathname,'/',files{1}];
    [finfo, acq, tran, ~, ~] = ReadRFE(filename, nr_load_fr, nr_skip_fr);
    frames = Frames(finfo,acq,tran,flow_settings,flow_rate,RF_settings,...
        pathname, objname);
    
    %load frames
    for i = 1:length(files);
        filename=[pathname,'/',files{i}];
        [~, ~, ~, ~, rfm] = ReadRFE(filename, nr_load_fr, nr_skip_fr);
        frames.ReadRFM(rfm);
    end
    
end


function [flow_settings, flow_rate, RF_settings] = LoadSettings(pathname)
currentpath=pwd;
cd(pathname)
msg='Choose Folder with getFlowSetting.m and getRFsettings.m';
while ~exist('getFlowSettings.m', 'file')
    disp(msg);
    disp(['Pathname = ',pathname]);
    temppath=uigetdir([pathname,'/..'],msg);
    cd(temppath);
end
[flow_settings, flow_rate] = getFlowSettings();

msg='Choose Folder with getRFsettings.m';
while ~exist('getRFsettings.m', 'file')
    disp(msg);
    temppath=uigetdir([pathname,'/..'],msg);
    cd(temppath);
end
RF_settings = getRFsettings();
cd(currentpath);
end