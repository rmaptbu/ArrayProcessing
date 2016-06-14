classdef Frames < handle
    properties
        %settings
        finfo %File info
        acq %Acquisition settings
        tran %Transmission settings     
        flw %Flow settings
        flw_r %Flow rate
        RF %RF settings
        pathname %location of data on disk
        filename %name of object when saved
        
        %data
        rfm %rf matrix
        p0_recon_TR
        p0_recon_FT
        
        %k-Wave settings
        kgrid
        medium
        source
        sensor
        input_args
        dt
        
        %Plots
        X %Ticks for X-axis scaling
        Y %Ticks for Y-axis scaling

    end
    methods
        function obj = Frames (finfo, acq, tran, flow_settings, ...
                flow_rate, RF_settings, pathname, filename)
            obj.finfo = finfo;
            obj.acq = acq;
            obj.tran = tran;
            obj.flw = flow_settings;
            obj.flw_r = flow_rate;
            obj.RF = RF_settings;
            obj.pathname = pathname; %location of data
            obj.filename = filename
            
            dx = obj.RF.pitch*1E-3;
            dy = obj.RF.speed_of_sound*obj.dt;
            obj.X = [0:dx:(obj.finfo.nrl-1)*dx]*1E3-(obj.finfo.nrl/2)*dx*1E3;
            obj.Y = [0:dy:(obj.finfo.nrs-1)*dy]*1E3;
        end
        function ReadRFM (obj, rfm)
            if obj.rfm
                obj.rfm = cat(3,obj.rfm,rfm);
            else
                obj.rfm = rfm;
            end
        end %Read RFM files
        function KWaveInit(obj) %Setup for k-wave toolbox
            c = obj.RF.speed_of_sound;
            obj.dt=1/obj.acq.fs;
            
            %nrl=number of lines, nrs=number of samples
            Nx = obj.finfo.nrl;  % number of grid points in the x (row) direction
            Ny = obj.finfo.nrs;  % number of grid points in the y (column) direction
            dx = obj.RF.pitch*1E-3;    % grid point spacing in the x direction [m]
            dy = c*obj.dt;         % grid point spacing in the y direction [m]
            obj.kgrid = makeGrid(Nx, dx, Ny, dy);
          
            obj.medium.sound_speed = c;	% [m/s]
            
            %create the time array
            obj.kgrid.t_array=0:obj.dt:((Ny-1)*obj.dt);
            
            %define a line sensor
            obj.sensor.mask = zeros(Nx, Ny);
            obj.sensor.mask(:, 1) = 1;
            
            %define the initial pressure
            obj.source.p0 = 0;
            
            PML_size = 20;
            obj.input_args = {'PMLInside', false, 'PMLSize', PML_size,...
                'Smooth', false};
            

        end
        function TR(obj,N) %Reconstruction via time reversal
            if ~N
                N=size(obj.rfm,3);
            end
            obj.p0_recon_TR=zeros(obj.finfo.nrs,obj.finfo.nrl,N);
            h = waitbar(0, 'Initialising Waitbar');
            msg='Computing time reversal...';
            for i=1:N
                waitbar(i/N,h,msg);
                sensor_data=obj.rfm(:,:,i)';
                obj.sensor.time_reversal_boundary_data = sensor_data;
                recon = kspaceFirstOrder2D(obj.kgrid, obj.medium,...
                    obj.source, obj.sensor, obj.input_args{:});
                % add first order compensation for only recording over a half plane
                recon = recon*2;
                obj.p0_recon_TR(:,:,i) = recon';
            end
            close(h);
            disp('Saving..');
            obj.save()
            disp('Done.');
        end
        function FT(obj,N) %Reconstruction via fourier transform
            if ~N
                N=size(obj.rfm,3);
            end
            dy = obj.RF.speed_of_sound*obj.dt;
            h = waitbar(0, 'Initialising Waitbar');
            msg='Computing FFT...';
            for i=1:N
                waitbar(i/N,h,msg);
                disp(i/N)
                sensor_data=obj.rfm(:,:,i);
                obj.p0_recon_FT(:,:,i) = kspaceLineRecon(sensor_data, dy, obj.dt, ...
                    obj.medium.sound_speed);%, 'Interp', '*linear');
            end
            close(h);
            disp('Saving..');
            obj.save()
            disp('Done.');
        end
        function save(obj) %save the object to original folder
            save([obj.pathname,'/',obj.filename,'.mat'],'obj')
        end
        function PlotRFM (obj)
            
            Im1=obj.rfm(10:end,:,1);
            Im2=obj.p0_recon_TR(10:end,:,1);
            Im3=obj.p0_recon_FT(10:end,:,1);
            figure;
            colormap('gray');
            
            subplot(1,3,1)
            imagesc(obj.X,obj.Y,Im1);
            title('Raw');
            caxis([-20 20])
            
            subplot(1,3,2)
            imagesc(obj.X,obj.Y,Im2);
            title('Time Reversal');
            caxis([-20 20])
            
            subplot(1,3,3)
            imagesc(obj.X,obj.Y,Im3);
            title('Fourier Transform');
            caxis([-20 20])
        end
        function PlotTR (obj) 
            Im1=obj.p0_recon_TR(100:end,:,1);
            Im2=obj.p0_recon_TR(100:end,:,2);
            figure;
            colormap('gray');
            
            subplot(1,2,1)
            imagesc(obj.X,obj.Y,Im1);
            title('Im1');
            
            subplot(1,2,2)
            imagesc(obj.X,obj.Y,Im2);
            title('Im2'); 
        end
    end
end

