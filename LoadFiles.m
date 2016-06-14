function [PAFrames, USFrames, pathname] = LoadFiles ()
%function settings
nr_load_fr=2; %how many consectives frames in file
nr_skip_fr=18; %how many frames to skip reading
%select folder
pathname = uigetdir('/Users/Thore/Documents/Transducer Measurements/160707LinArray');

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
filename=[pathname,'/',PAFiles{1}];
[finfo, acq, tran, ~, ~] = ReadRFE(filename, nr_load_fr, nr_skip_fr);
PAFrames = Frames (finfo,acq,tran,flow_settings,flow_rate,RF_settings);

%load frames
for i = 1:length(PAFiles);
    filename=[pathname,'/',PAFiles{i}];
    [~, ~, ~, ~, rfm] = ReadRFE(filename, nr_load_fr, nr_skip_fr);
    PAFrames.ReadRFM(rfm);
end

%% Load US data
% load file settings
filename=[pathname,'/',USFiles{1}];
[finfo, acq, tran, ~, ~] = ReadRFE(filename, 20, 0);
USFrames = Frames (finfo,acq,tran,flow_settings,flow_rate,RF_settings);

%load frames
for i = 1:length(USFiles);
    filename=[pathname,'/',USFiles{i}];
    [~, ~, ~, ~, rfm] = ReadRFE(filename, 20, 0);
    USFrames.ReadRFM(rfm);
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