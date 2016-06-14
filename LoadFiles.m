function [PAFrames, USFrames, pathname] = LoadFiles(varargin)
%function settings
nr_load_fr=2; %how many consectives frames in file
nr_skip_fr=18; %how many frames to skip reading

if ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'Path'
                pathname = varargin{input_index + 1};
            otherwise
                error('Unknown optional input');
        end
    end
end              
%select folder
if ~exist('pathname','var')
    pathname = uigetdir('/Users/Thore/Documents/Transducer Measurements/160707LinArray');
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

%% Load PA data
% load file settings
PAFrames = LoadFrames(pathname, 'PAFrames', PAFiles, nr_load_fr, nr_skip_fr);

%% Load US data
% load file settings
filename=[pathname,'/',USFiles{1}];
[finfo, acq, tran, ~, ~] = ReadRFE(filename, 20, 0);
USFrames = Frames (finfo,acq,tran,flow_settings,flow_rate,RF_settings,...
    pathname, 'USFrames');

%load frames
for i = 1:length(USFiles);
    filename=[pathname,'/',USFiles{i}];
    [~, ~, ~, ~, rfm] = ReadRFE(filename, 20, 0);
    USFrames.ReadRFM(rfm);
end
PAFrames.KWaveInit();
USFrames.KWaveInit();
end

function [frames] = LoadFrames(pathname,objname,files,nr_load_fr, nr_skip_fr)
currentpath=pwd;
cd(pathname)
if exist([objname,'.mat'], 'file')
    title=['Load ',objname];
    msg=['Load existing ',objname,'.mat?'];
    button=questdlg(msg,title,'Yes','No','Yes');
    switch button
        case 'Yes'
            frames = load([pathname,'/',objname,'.mat']);
            frames = frames.obj;
            cd(currentpath)
            return
        case 'No'
            % load file settings
    end
end
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
    cd(currentpath)
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