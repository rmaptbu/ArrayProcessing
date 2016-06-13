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
        
        %k-Wave data
        kgrid
        medium
        source
        sensor
        input_args
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
        function kWav(obj) %Setup for k-wave toolbox
            PML_size = 20;          % size of the PML in grid points
            Nx = 128 - 2*PML_size;  % number of grid points in the x (row) direction
            Ny = 256 - 2*PML_size;  % number of grid points in the y (column) direction
            dx = 0.1e-3;            % grid point spacing in the x direction [m]
            dy = 0.1e-3;            % grid point spacing in the y direction [m]
            obj.kgrid = makeGrid(Nx, dx, Ny, dy);
            
            obj.medium.sound_speed = obj.RF.speed_of_sound;	% [m/s]

        end
        function TR(obj) %Reconstruction through time reversal
            p0_recon = kspaceFirstOrder2D(kgrid, medium, source, sensor, input_args{:});
            
        end
        function PlotRFM (obj)
            Im_mean=mean(obj.rfm,3);
            imagesc(Im_mean);
        end
    end
end

