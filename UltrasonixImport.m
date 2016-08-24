clear all
close all
addpath('D:\Wenfeng\k_wave_toolbox_version_1.0\k_Wave_Toolbox');

addpath('D:\SonixProgramming\sdk610\MATLAB\SonixDataTools\DataReaders\DAQ')
% addpath('O:\Original\sdk610\MATLAB\SonixDataTools\GUI')
addpath('D:\SonixProgramming\sdk610\MATLAB\SonixDataTools\DataReaders\RPread')
% Load data

path1 = 'D:\DAQRT\1\';
% path2 = 'D:\DAQRT\2\';

cd('D:\Wenfeng\sdk610\stage\bin\Release')
[Im,header] = RPread('TexoData.rf');


aver = mean(Im, 3);
cd('D:\DAQRT\')
%
% [Im,header] = RPread('TexoData_transmitting.rf');


samplef = 40e6;

finerecon = 1;
% create the computational grid
PML_size = 20;          % size of the PML in grid points
Ny = 128*finerecon;  % number of grid points in the x (row) direction
Nx = 256*finerecon;  % number of grid points in the y (column) direction
dx = 0.305e-3./finerecon;            % grid point spacing in the x direction [m]
dy = 0.305e-3./finerecon;            % grid point spacing in the y direction [m]
kgrid = makeGrid(Nx, dx, Ny, dy);

% define the properties of the propagation medium
medium.sound_speed = 1490;           % [m/s]

%  T=22;
%  medium.sound_speed=1.40238744e3+5.03836171*T-5.81172916e-2*(T.^2)+3.34638117e-4*(T^3)-...
%      1.48259672e-6*(T^4)+3.16585020e-9*(T^5);

% create initial pressure distribution using makeDisc

% define a binary line sensor
sensor.mask = zeros(Nx, Ny);
sensor.mask(1, 1:finerecon:Ny) = 1;

% create the time array
dt = 1/samplef;

% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
input_args = {'PMLInside', false, 'PMLSize', PML_size, 'PlotPML', false, 'Smooth', false};


[header, data] = readDAQ(path1, ones(1,128), 2, true);
kgrid.t_array = (1:max(size(data))).*dt;
frameNo = header(2);
Recon_data = zeros(header(3)+160,(header(1)+1)*finerecon,header(2));
data1 = zeros(header(3)+160,(header(1)+1));
for i=1:frameNo;
    [header, data] = readDAQ(path1, ones(1,128), i, true);
    data1(161:end,:) =  data;
    data1(1:200,:) = 0;
    
    sensor_data = data1.';
    X = 1:finerecon:Ny;
    X1 = 1:1:Ny;
    Y1 = 1:max(size(sensor_data));
    Y1 = Y1';
    sensor_data2 = interp2(X,Y1,sensor_data.',X1,Y1,'nearest');
    sensor_data2(:,end)=0;
    %sensor_data2 = sensor_data2.';
    
    % reconstruct the initial pressure
    % FFT recon
    p_xy = kspaceLineRecon(sensor_data2, dy, dt, medium.sound_speed, 'Plot', false, 'PosCond', true);
    %  p = recon_DelayAndSum_line(sensor_data2, dy, dt, medium.sound_speed, 'Plot', true);
    
    
    Recon_data(:,:,i) = p_xy;
    
    % time-reversal
    %  source.p0 = 0;
    %  sensor.time_reversal_boundary_data = sensor_data./max(max(sensor_data));
    % p0_recon = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
    %inten_aver(i,k) =  mean(mean(abs(Recon_data(1170:1250,250:290,i))));  % blood
    % inten_aver(i,k) =  mean(mean(abs(Recon_data(800:840,160:230,i)))); % fat
end

%  p_xy(p_xy < 0) = 0;
p_xy = squeeze(mean(Recon_data(:,:,:),3));
%p_xy = squeeze((Recon_data(:,:,1)));

p_xy = abs(hilbert(squeeze(p_xy)));
datamax = max(p_xy(:));
p_xy_dB = 20*log10(p_xy./datamax);
%   p_xy = p_xy./max(max(p_xy));





figure('units','normalized','outerposition',[0 0 1 1])
h1=subplot(2,2,1);
imagesc([1:Ny]*dy*1e3, [1:max(size(sensor_data2))]*dt*medium.sound_speed*1e3, p_xy_dB(1:end,:));

set(h1,'FontSize',20, 'TickDir','out','TickLength',[.01 0]);
axis image;
xlabel(['X [mm]'],'FontSize',20);
ylabel(['Z [mm]'],'FontSize',20);
xlim([0 38]);
ylim([0 70]);
colorbar('FontSize',20);
caxis([-45,0])
colormap(h1, hot);
colorbar('off')
title('PA','FontSize',20)



%   p_xy = abs(hilbert(Im(:,:,1)));
us_xy = 20*log10(abs(hilbert(aver))./max(max(abs(hilbert(aver)))));
%  p_xy = aver;
% h2=subplot(2,3,2);
% imagesc([3:258]*1.5e-3*100 - 0.6, [1:max(size(sensor_data2))]*dt*medium.sound_speed*1e3/2, us_xy(1:end,3:end));
% 
% set(h2,'FontSize',20, 'TickDir','out','TickLength',[.01 0]);
% axis image;
% xlabel(['X [mm]'],'FontSize',20);
% ylabel(['Z [mm]'],'FontSize',20);
% xlim([0 38]);
% ylim([0 40]);
% colorbar('FontSize',20);
% caxis([-70,0])
% colormap(h2, gray);

% colorbar('off')
% title('US','FontSize',20)



subplot(2,2,2)
us_xy_norm = (us_xy + 60)./ max(max(us_xy + 60));

subimage([1:Ny]*dy*1e3, [1:max(size(sensor_data2))]*1/80e6*medium.sound_speed*1e3, us_xy_norm(1:end,3:end));
set(gca,'FontSize',20, 'TickDir','out','TickLength',[.01 0]);
axis image;
xlabel(['X [mm]'],'FontSize',20);
ylabel(['Z [mm]'],'FontSize',20);
xlim([0 38]);
ylim([0 40]);
colorbar('FontSize',20);
% caxis([0,1])
colormap(gray);colorbar('off')


title('PA + US','FontSize',20)
hold on;

hImg = imagesc([1:Ny]*dy*1e3, [1:max(size(sensor_data2))]*dt*medium.sound_speed*1e3/2, p_xy_dB(1:end,:));
colormap hot;  colorbar('off')
caxis([-35,0])
set(hImg, 'AlphaData', 0.3);


for i=1:frameNo;
    [header, data] = readDAQ(path2, ones(1,128), i, true);
    data1(161:end,:) =  data;
    data1(1:250,:) = 0;
    
    sensor_data = data1.';
    X = 1:finerecon:Ny;
    X1 = 1:1:Ny;
    Y1 = 1:max(size(sensor_data));
    Y1 = Y1';
    sensor_data2 = interp2(X,Y1,sensor_data.',X1,Y1,'nearest');
    sensor_data2(:,end)=0;
    %sensor_data2 = sensor_data2.';
    
    % reconstruct the initial pressure
    % FFT recon
    p_xy = kspaceLineRecon(sensor_data2, dy, dt, medium.sound_speed, 'Plot', false, 'PosCond', true);
    %  p = recon_DelayAndSum_line(sensor_data2, dy, dt, medium.sound_speed, 'Plot', true);
    
    
    Recon_data(:,:,i) = p_xy;
    
    % time-reversal
    %  source.p0 = 0;
    %  sensor.time_reversal_boundary_data = sensor_data./max(max(sensor_data));
    % p0_recon = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
    %inten_aver(i,k) =  mean(mean(abs(Recon_data(1170:1250,250:290,i))));  % blood
    % inten_aver(i,k) =  mean(mean(abs(Recon_data(800:840,160:230,i)))); % fat
end

%  p_xy(p_xy < 0) = 0;
p_xy = squeeze(mean(Recon_data(:,:,:),3));
%p_xy = squeeze((Recon_data(:,:,1)));
p_xy = abs(hilbert(squeeze(p_xy)));
datamax = max(p_xy(:));
p_xy_dB = 20*log10(p_xy./datamax);
%   p_xy = p_xy./max(max(p_xy));


%figure('units','normalized','outerposition',[0 0 1 1])
h4=subplot(2,2,3);
imagesc([1:Ny]*dy*1e3, [1:max(size(sensor_data2))]*dt*medium.sound_speed*1e3, p_xy_dB(1:end,:));

set(h4,'FontSize',20, 'TickDir','out','TickLength',[.01 0]);
axis image;
xlabel(['X [mm]'],'FontSize',20);
ylabel(['Z [mm]'],'FontSize',20);
xlim([0 38]);
ylim([0 60]);
colorbar('FontSize',20);
caxis([-25,0])
colormap(h4, hot);
colorbar('off')
title('Wavelength: 1300 nm','FontSize',20)



% %   p_xy = abs(hilbert(Im(:,:,1)));
% us_xy = 20*log10(abs(hilbert(aver))./max(max(abs(hilbert(aver)))));
% %  p_xy = aver;
% h5=subplot(2,3,5);
% imagesc([3:258]*1.5e-3*100 - 0.6, [1:max(size(sensor_data2))]*dt*medium.sound_speed*1e3/2, us_xy(1:end,3:end));
% 
% set(h5,'FontSize',20, 'TickDir','out','TickLength',[.01 0]);
% axis image;
% xlabel(['X [mm]'],'FontSize',20);
% ylabel(['Z [mm]'],'FontSize',20);
% xlim([0 38]);
% ylim([0 40]);
% colorbar('FontSize',20);
% caxis([-70,0])
% colormap(h5, gray);
% 
% colorbar('off')
% title('US','FontSize',20)



subplot(2,2,4)
us_xy_norm = (us_xy + 60)./ max(max(us_xy + 60));

subimage([1:Ny]*dy*1e3, [1:max(size(sensor_data2))]*dt*medium.sound_speed*1e3/2, us_xy_norm(1:end,3:end));
set(gca,'FontSize',20, 'TickDir','out','TickLength',[.01 0]);
axis image;
xlabel(['X [mm]'],'FontSize',20);
ylabel(['Z [mm]'],'FontSize',20);
xlim([0 38]);
ylim([0 40]);
colorbar('FontSize',20);
% caxis([0,1])
colormap(gray);colorbar('off')


title('PA + US','FontSize',20)
hold on;

hImg = imagesc([1:Ny]*dy*1e3, [1:max(size(sensor_data2))]*dt*medium.sound_speed*1e3, p_xy_dB(1:end,:));
colormap hot;  colorbar('off')
caxis([-35,0])
set(hImg, 'AlphaData', 0.3);