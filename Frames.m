classdef Frames < handle
    properties
        %settings
        finfo %File info
        acq %Acquisition settings
        tran %Transmission settings     
        flw %Flow settings
        flw_r %Flow rate
        RF %RF settings
        
        %data
        rfm %rf matrix
        p0_recon_TR
        p_xy
        
        %k-Wave settings
        kgrid
        medium
        source
        sensor
        input_args
        dt

    end
    methods
        function obj = Frames (finfo, acq, tran, flow_settings, ...
                flow_rate, RF_settings)
            obj.finfo = finfo;
            obj.acq = acq;
            obj.tran = tran;
            obj.flw = flow_settings;
            obj.flw_r = flow_rate;
            obj.RF = RF_settings;
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
            obj.sensor.mask(1, :) = 1;
            
            %define the initial pressure
            obj.source.p0 = 0;
            
            PML_size = 20;
            obj.input_args = {'PMLInside', false, 'PMLSize', PML_size,...
                'Smooth', false, 'PlotPML', false};
            

        end
        function TR(obj) %Reconstruction through time reversal

            sensor_data=obj.rfm(:,:,1);
            obj.sensor.time_reversal_boundary_data = sensor_data;
            obj.p0_recon_TR = kspaceFirstOrder2D(obj.kgrid, obj.medium,...
                obj.source, obj.sensor, obj.input_args{:});
            % add first order compensation for only recording over a half plane
            obj.p0_recon_TR = obj.p0_recon_TR*2;
            
            % repeat the FFT reconstruction for comparison
%             obj.p_xy = kspaceLineRecon(sensor_data.', obj.dy, obj.dt, ...
%                 obj.medium.sound_speed, 'PosCond', true, 'Interp', '*linear');
            
        end
        function PlotRFM (obj)
            Im_mean=mean(obj.rfm,3);
            imagesc(Im_mean);
        end
    end
end

