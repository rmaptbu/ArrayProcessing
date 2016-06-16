classdef Frames < handle
    %Object to handle 2D Photoacoustic Frames
    %Thore Bucking 2016
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
        rfm_b %backup of rfm matrix
        p0_recon_TR %reconstruction using time reversal
        p0_recon_FT %reconstruction using fourier transform
        QS1 %number of initial points to discard: Frame 1
        QS2 %number of initial points to discard: Frame 2
        
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
            obj.filename = filename; %name of .mat file
            
            dx = obj.RF.pitch*1E-3;
            dy = obj.RF.speed_of_sound/obj.acq.fs;
            obj.X = [0:dx:(obj.finfo.nrl-1)*dx]*1E3-(obj.finfo.nrl/2)*dx*1E3;
            obj.Y = [0:dy:(obj.finfo.nrs-1)*dy]*1E3;
            
            obj.QS1=0;
            obj.QS2=0;
        end
        function ReadRFM (obj, rfm)
            if obj.rfm
                obj.rfm = cat(3,obj.rfm,rfm);
            else
                obj.rfm = rfm;
            end
            obj.rfm_b = obj.rfm; %backup original rf matrix
        end %Read RFM files
        function LoadRFM(obj) %
            if ~isempty(obj.rfm_b)
                obj.rfm=obj.rfm_b;
            else
                warning('There is no backup RF data stored.')
            end
        end %Load RF data from internal backup
        function KWaveInit(obj, varargin) %Setup for k-wave toolbox
            upsample=1;
            if ~isempty(varargin)
                for input_index = 1:2:length(varargin)
                    switch varargin{input_index}
                        case 'Upsample'
                            upsample = varargin{input_index + 1};
                        otherwise
                            error('Unknown optional input');
                    end
                end
            end
            
            
            c = obj.RF.speed_of_sound;
            obj.dt=1/obj.acq.fs;
            
            %nrl=number of lines, nrs=number of samples
            Nx = obj.finfo.nrl;  % number of grid points in the x (row) direction
            Ny = size(obj.rfm,1);  % number of grid points in the y (column) direction
            dx = obj.RF.pitch*1E-3;    % grid point spacing in the x direction [m]
            dy = c*obj.dt;         % grid point spacing in the y direction [m]
            
            if upsample
            %upsampling detectors to match dy...
            N = ceil(dx/dy);
            dx = dx/N;
            Nx = Nx*N;
            obj.Upsample(N);
            end
            
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
            obj.p0_recon_TR=zeros(size(obj.rfm,1),size(obj.rfm,2),N);
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
            obj.Save()
            disp('Done.');
        end
        function FT(obj,N,varargin) %Reconstruction via fourier transform
            if ~N
                N=size(obj.rfm,3);
            end
            paddingX=size(obj.rfm,2); %zero padding outside of frame to prevent wrapping
            paddingY=0;
            save_opt = 0;
            if ~isempty(varargin)
                for input_index = 1:2:length(varargin)
                    switch varargin{input_index}
                        case 'Save'
                            save_opt = varargin{input_index + 1};
                        case 'Padding'
                            paddingY = varargin{input_index + 1};
                        otherwise
                            error('Unknown optional input');
                    end
                end
            end
            dy = obj.RF.speed_of_sound*obj.dt;
            obj.p0_recon_FT=zeros(size(obj.rfm,1),size(obj.rfm,2),N);
            h = waitbar(0, 'Initialising Waitbar');
            msg='Computing FFT...';
            dimX=size(obj.rfm(:,:,1),2);
            dimY=size(obj.rfm(:,:,1),1);
            
            sensor_data=zeros(2*paddingY+dimY,2*paddingX+dimX);
            for i=1:N
                waitbar(i/N,h,msg);
                disp(i/N)     
                sensor_data(paddingY+1:paddingY+dimY,paddingX+1:paddingX+dimX)=obj.rfm(:,:,i);
                recon = kspaceLineRecon(sensor_data, dy, obj.dt, ...
                    obj.medium.sound_speed, 'Interp', '*linear');  
                
                obj.p0_recon_FT(:,:,i) = recon(paddingY+1:paddingY+dimY,paddingX+1:paddingX+dimX);
            end
            close(h);
            if save_opt
                disp('Saving..');
                obj.Save()
                disp('Done.');
            end
        end
        function [QS1, QS2] = QSCorrect(obj,QS1, varargin) %Q Switch correction
            %Removes initial data points of frame pairs
            %Only necessary if acquisition is not triggered by the laser
            %output
            if isempty(QS1) %accept empty input to atomatically detect delay
                disp('Auto detect QS delay');
                %find first point where there is no background
                %thresholded by n sigma
                n=5;
                flag1=0;
                flag2=0;
                for i=5:size(obj.rfm_b,1)
                    mean_fr1=mean(obj.rfm_b(1:i,end,1:2:end-1),3);
                    mean_fr2=mean(obj.rfm_b(1:i,end,2:2:end),3);
                    sig1=std(mean_fr1(1:end-1));
                    sig2=std(mean_fr2(1:end-1));
                    
                    thresh1=n*sig1+abs(mean_fr1(end-1));
                    thresh2=n*sig2+abs(mean_fr2(end-1));
                    if abs(mean_fr1(end))>thresh1 && ~flag1
                        flag1 = 1;
                        obj.QS1=i;
                        QS1=i/(1E-9*obj.acq.fs);
                    end
                    if abs(mean_fr2(end))>thresh2 && ~flag2
                        flag2 = 1;
                        obj.QS2=i;
                        QS2=i/(1E-9*obj.acq.fs);
                    end
                    if flag1 && flag2
                        break
                    end
                end
            else   
                %QS1 and QS2 given in ns
                %Convert ns to data points to discard
                obj.QS1=QS1*1E-9*obj.acq.fs;
                if ~isempty(varargin)
                    QS2=QS1;
                else                    
                    QS2=varargin{1};
                end
                obj.QS2=QS2*1E-9*obj.acq.fs;
            end
            
            h = waitbar(0, 'Initialising Waitbar');
            msg='Resampling...';
            obj.rfm=[];
            for i=1:2:size(obj.rfm_b,3)-1
            waitbar(i/(size(obj.rfm_b,3)-1),h,msg);
            obj.rfm(:,:,i)=imtranslate(obj.rfm_b(:,:,i),[0 -obj.QS1]);
            obj.rfm(:,:,i+1)=imtranslate(obj.rfm_b(:,:,i+1),[0 -obj.QS2]);
            end
            obj.rfm=obj.rfm(1:end-floor(obj.QS2),:,:);
            
            %update Y tick labels
            dy = obj.RF.speed_of_sound*obj.dt;
            obj.Y = [0:dy:(size(obj.rfm,1)-1)*dy]*1E3;
            
            close(h);      
            
        end
        function RemoveNoise(obj,L) %input in mm
            %Number of initial points to set to 0
            dy = obj.RF.speed_of_sound/obj.acq.fs*1E3;
            N = floor(L/dy);
            obj.rfm(1:N,:,:)=0;
        end
        function Upsample(obj,N) %Upsample number of detectors
            rfm_t = obj.rfm;
            for i=1:size(rfm_t,2)
                for j=1:N
                    %nearest neighbour interpolation
                    obj.rfm(:,N*(i-1)+j,:)=rfm_t(:,i,:);
                end
            end
        end
        function Save(obj) %save the object to original folder
            save([obj.pathname,'/',obj.filename,'Del',num2str(obj.QS1),'.mat'],'obj')
        end
        function PlotRFM (obj, varargin)
            
            SaveFig = 0;
            Average = 0;
            figname = 0;
            if ~isempty(varargin)
                for input_index = 1:2:length(varargin)
                    switch varargin{input_index}
                        case 'SaveFig'
                            SaveFig = varargin{input_index + 1};
                        case 'FigName'
                            figname = varargin{input_index + 1};
                        case 'Average'
                            Average = varargin{input_index + 1};
                        otherwise
                            error('Unknown optional input');
                    end
                end
            end
            
            if Average
                Im1=mean(obj.rfm,3);
                Im2=mean(obj.p0_recon_TR,3);
                Im3=mean(obj.p0_recon_FT,3);
            else
                Im1=obj.rfm(:,:,1);
                Im2=obj.p0_recon_TR(:,:,1);
                Im3=obj.p0_recon_FT(:,:,1);
            end
            
            fig=figure;
            colormap('gray');
            
            subplot(1,3,1)
            imagesc(obj.X,obj.Y,Im1);
            title('Raw');
            caxis([-40 40])
            xlabel('Lateral (mm)');
            ylabel('Depth (mm)');
            
            subplot(1,3,2)
            imagesc(obj.X,obj.Y,Im2);
            title('Time Reversal');
            caxis([-80 80])
            xlabel('Lateral (mm)');
            ylabel('Depth (mm)');            
            
            subplot(1,3,3)
            imagesc(obj.X,obj.Y,Im3);
            title('Fourier Transform');
            caxis([-40 40])
            xlabel('Lateral (mm)');
            ylabel('Depth (mm)');
            
            set(fig, 'Position', [100 100 800 600]);
            
            if SaveFig
                if ~figname
                    figname = [obj.pathname,'/',obj.filename,'Del',num2str(obj.QS1),'.png'];
                else
                    figname = [obj.pathname,'/',figname,'.png'];
                end
                set(gcf,'PaperPositionMode','auto')
                print(fig,figname,'-dpng','-r0')
                close(fig);
            end
        end
        function PlotTR (obj) 
            Im1=obj.p0_recon_TR(:,:,1);
            Im2=obj.p0_recon_TR(:,:,2);
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

