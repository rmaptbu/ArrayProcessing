function reconObj = ReconFcn(rfObj,varargin)
%% Parse Inpus
upsample=2;
recon_type='FT';
if ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'Upsample'
                upsample = varargin{input_index + 1};
            case 'ReconType'
                recon_type = varargin{input_index + 1};
            otherwise
                error('Unknown optional input');
        end
    end
end

%% Get RF Matrix
c = rfObj.c; %m/s
dt=1/rfObj.fs; %s
dx = rfObj.pitch;    % grid point spacing in the x direction [m]
dy = c*dt;         % grid point spacing in the y direction [m]

del=rfObj.del; %in mm
del=del/(dy*1E3); %in pixel
h = waitbar(0, 'Initialising Waitbar');
msg='Extracting RF Matrix...';
rfm=rfObj.rfm;
for i=1:size(rfm,3)
    waitbar(i/(size(rfObj,3)),h,msg);
    rfm(:,:,i)=imtranslate(rfObj.rfm(:,:,i),[0 -del]);  
end
rfm=rfm(1:end-floor(del),:,:);

%update Y tick labels
Y = (0:dy:(size(rfm,1)-1)*dy)*1E3;
X = rfObj.X;
close(h);


%% Initialise KWave
Nx = size(rfm,2);  % number of grid points in the x (row) direction
Ny = size(rfm,1);  % number of grid points in the y (column) direction
Nz = size(rfm,3);
%Upsampling
dx = dx/upsample;
[x,y]=meshgrid(X,Y);
X=interp1(1:length(X),X,1:1/upsample:length(X));
[xq,yq]=meshgrid(X,Y);
Nx = size(xq,2);
rfm_t=zeros(Ny,Nx,Nz);
for i=1:size(rfm,3)
    rfm_t(:,:,i)=interp2(x,y,rfm(:,:,i),xq,yq,'nearest');
end
rfm=rfm_t;
clear('x','y','xq','yq','rfm_t');

kgrid = makeGrid(Nx, dx, Ny, dy);
medium.sound_speed = c;	% [m/s]
%create the time array
kgrid.t_array=0:dt:((Ny-1)*dt);
%define a line sensor
sensor.mask = zeros(Nx, Ny);
sensor.mask(:, 1) = 1;
%define the initial pressure
source.p0 = 0;

PML_size = 20; %Perfectly matched layer
input_args = {'PMLInside', false, 'PMLSize', PML_size,...
    'Smooth', false, 'PlotSim',true};

%% Reconstruction
switch(recon_type)
    case 'TR'
        TR();
    case 'FT'
        FT();
    case 'BF'
        BF();
    otherwise
        error('Unknown optional input');
end
%% Create ouptut object
reconObj = ReconObj(rfObj.fs, rfObj.c, rfObj.pitch,...
    rfObj.pathname,'recon', X, Y,p0_recon);

%% Reconstruction Functions
    function TR() %Reconstruction via time reversal
        p0_recon=zeros(Ny,Nx,Nz);
        h = waitbar(0, 'Initialising Waitbar');
        msg='Computing time reversal...';
        
        for j=1:Nz
            waitbar(j/Nz,h,msg);
            sensor_data=rfm(:,:,j)';
            sensor.time_reversal_boundary_data = sensor_data;
            recon = kspaceFirstOrder2D(kgrid, medium,...
                source, sensor,input_args{:});
            % add first order compensation for only recording over a half plane
            recon = recon*2;
           p0_recon(:,:,j) = recon';
        end
        close(h);

    end
    function FT() %Reconstruction via fourier transform
        pad=floor(1.5*Nx); %zero padding outside of frame to prevent wrapping        
        p0_recon=zeros(Ny,Nx,Nz);
        h = waitbar(0, 'Initialising Waitbar');
        msg='Computing FFT...';        
        sensor_data=zeros(Ny,2*pad+Nx); 
        
        for j=1:Nz
            waitbar(j/Nz,h,msg);
%             disp(j/Nz);
            sensor_data(1:Ny,pad+1:pad+Nx)=rfm(:,:,j);
            recon = kspaceLineRecon(sensor_data, dx, dt, ...
               c, 'Interp', '*linear');
            
            p0_recon(:,:,j) = recon(1:Ny,pad+1:pad+Nx);
        end
        close(h);
    end
    function BF(varargin) %Reconstruction via beam forming
        %Code adapted from Pim van den Berg (University of Twente)
        if ~isempty(varargin)
            aptr=varargin{1};
        else
            aptr=Nx;
        end
        aptr=ceil((aptr-1)/2.0);
        
        p0_recon=zeros(Ny,Nx,Nz);

        pitch = dx.* rfObj.fs/c;
        
        h = waitbar(0, 'Initialising Waitbar');
        msg='Computing Beam Forming...';
        for j=1:Nz
            data=squeeze(rfm(:,:,j));
            waitbar(j/Nz,h,msg);
            for channel_i = 1:Nx % for every channel
                result_xx = zeros(1,Ny);
                for d_y = 1:Ny % for all rows
                    for d_x = 1:Nx
                        if(dx ~= channel_i) % for all off_axis image points
                            if abs(channel_i - d_x) <= aptr
                                delta_x = (channel_i - d_x) * pitch; % Pythagoras
                                delta_y = d_y;
                                % distance from imagepoint to element:
                                l = sqrt(delta_x^2 + delta_y^2);
                                l = round(l); % rounding distance
                                % to whole sample point,
                                % no interpolation
                                % check that distance is within recon area:
                                if(l >= 1 && l <= Ny)
                                    result_xx(d_y) = result_xx(d_y) + data(l,d_x); % sum image points
                                end
                            end
                        else % for the on-axis image points
                            result_xx(d_y) = result_xx(d_y) + data(d_y, channel_i); % sum image points
                        end
                    end
                end
                p0_recon(:, channel_i, j) = result_xx;
            end
        end
        close(h);

        
    end

end

