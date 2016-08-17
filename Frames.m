classdef Frames < handle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Object to handle 2D Photoacoustic Frames
    %Thore Buecking 2016 (rmaptbu@ucl.ac.uk)
    %Requires the K-Wave Toolbox (www.k-wave.org)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Example usage:
    %1.)obj.ReadRFM(rfm) >> Read Data rfm(rows,col,frames)
    %   |obj.QSCorrect([]) >> Autodetect Laser firing
    %2.)|OR (if timing is already correct)
    %   |obj.KWaveInit() >> Prepare K-Wave settings 
    %3.)obj.FT(0) >> reconstruct all frames using fourier transfrom  
    %4.)obj.PlotRFM() >> Make sure Reconstruction worked
    %5.)obj.EnsembleCorrelation('FT') >> Create Cross-Correlation Map
    %6.)obj.PlotXC() >> View Result
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Methods
    %Frames: Initialise Object
    %       -RF, finfo, acq, tran: Files from Data acquisition
    %           uses: obj.RF.speed_of_sound
    %                 obj.RF.pitch: Pitch between detectors (mm)
    %                 obj.finfo.nrl: Number of lines (i.e. detectors)                    
    %                 obj.acq.fs: Sambling Frequency (Hz)
    %ReadRFM: Load RF matrix (rfm). Attaches to rfm and copies to rfm_b.
    %LoadRFM: Override rfm with rfm_b.
    %KWaveInit('Upsample',N): Create K-wave settings. Upsamples number of
    %       detectors in rfm. If N is not specified, it chooses N to be
    %       closest to making the grid square.
    %TR(N): Time reversal reconstruction. N=number of frames to compute.
    %       N=0: compute all frames.
    %FT(N, 'Padding', P): Fourier transform reconstruction. N=number of 
    %       frames to compute. N=0: compute all frames. P is number of
    %       pixels padded outside of frame to prevent wrapping artefacts.
    %       Default values should be good enough.
    %QSCorrect(QS1,QS2): Remove first QS(ns) from acquisition frame.
    %       If QS2 is not specified, QS2=QS1.
    %       If QS1=[], laser shot will be automatically detected
    %       !!!>>>Overrides rfm from rfm_b<<<!!! It will do KWaveInit()!
    %RemoveNoise(L): Replace first L(mm) points of data with zeros
    %Upsample(N): Upsamples number of detectes(multiplies by N)
    %EnsembleCorrelation(obj,type): type='FT' or 'TR'.Ensemble Correlate 
    %       all frames of type (TR reconstr. or FT reconstr.). Saves output 
    %       in xc_raw 
    %FindShift(): Take xc_raw and extract displacement and amplitude
    %       information in xc_disp and xc_amp (displacement and amplitude 
    %       respectively).
    %Save(): Save object to folder specified in the construction.
    %PlotRFM('SaveFig','FigName','Average'): Plot data, TR, FT.
    %       'SaveFig':Close figure, save in original folder
    %       'FigName':Specify name to save figure
    %       'Average':Plot average of all frames instead of first frame
    %PlotXC('SaveFig','FigName'): Plot Reconstrct, Xcorr_disp, Xcorr_amp.
    %       'SaveFig':Close figure, save in original folder
    %       'FigName':Specify name to save figure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        p0_recon
        p0_recon_TR %reconstruction using time reversal
        p0_recon_FT %reconstruction using fourier transform
        p0_recon_BF %reconstruction using beam forming
        QS1 %number of initial points to discard: Frame 1
        QS2 %number of initial points to discard: Frame 2
        
        %cross correlation results
        xc_raw %Raw ensemble cross correlation
        xc_raw_min
        xc_disp %cross correlation displacement map
        xc_amp %cross correlation amplitude map
        xc_amp_me
        xc_mask
        xc_flw %calculated flow speed
        xc_flw_std %standard deviation of flow speed
        xc_all %all flow speed measurements
        
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
        Y_xc %Ticks for Y-axis scaling of xcorrs
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
            obj.X = (0:dx:(obj.finfo.nrl-1)*dx)*1E3-(obj.finfo.nrl/2)*dx*1E3;
            obj.Y = (0:dy:(obj.finfo.nrs-1)*dy)*1E3;
            
            obj.QS1=0;
            obj.QS2=0;
        end
        %Prepare Object
        function ReadRFM (obj, rfm)
            if ~isempty(obj.rfm)
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
            c = obj.RF.speed_of_sound;
            obj.dt=1/obj.acq.fs;
            
            %nrl=number of lines, nrs=number of samples
            Nx = size(obj.rfm,2);  % number of grid points in the x (row) direction
            Ny = size(obj.rfm,1);  % number of grid points in the y (column) direction
            dx = obj.RF.pitch*1E-3;    % grid point spacing in the x direction [m]
            dx = dx/(Nx/size(obj.rfm_b,2));  % account for previous upsampling
            dy = c*obj.dt;         % grid point spacing in the y direction [m]
                      
            q=1;
            if ~isempty(varargin)
                if isnumeric(varargin{1});
                    Nx = varargin{1};
                    Ny = varargin{2};                
                    q=3;
                end
                for input_index = q:2:length(varargin)
                    switch varargin{input_index}
                        case 'Upsample'
                            upsample = varargin{input_index + 1};
                        otherwise
                            error('Unknown optional input');
                    end
                end
            end
                       
            %Upsampling
            %If not specified, choose value to make grid as close to square
            %as possible
            if exist('upsample','var')
                N = upsample;
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
            
            PML_size = 20; %Perfectly matched layer
            obj.input_args = {'PMLInside', false, 'PMLSize', PML_size,...
                'Smooth', false};
            

        end
        function [QS1, QS2] = QSCorrect(obj,varargin) %Q Switch correction
            %Removes initial data points of frame pairs
            %Only necessary if acquisition is not triggered by the laser
            %output
            auto = 0;
            if ~isempty(varargin)
                for input_index = 1:2:length(varargin)
                    switch varargin{input_index}
                        case 'AutoDetect'
                            auto = varargin{input_index + 1};
                        case 'QS1'
                            assert(~auto);
                            QS1 = varargin{input_index + 1}; 
                            obj.QS1=QS1*1E-9*obj.acq.fs;
                        case 'QS2'
                            assert(~auto);
                            QS2 = varargin{input_index + 1};
                            obj.QS2=QS2*1E-9*obj.acq.fs;
                        otherwise
                            error('Unknown optional input');
                    end
                end
            else
                TrigDel = obj.RF.laser_trig_delay + obj.RF.frame_skip + obj.RF.client_delay;
                AcqDel = obj.acq.dt/obj.acq.fs*1e6; %?s
                QS1 = (AcqDel - TrigDel)*obj.acq.fs/1e6;
                QS2 = QS1;
                obj.QS2 = QS2;
                obj.QS1 = QS1;
            end
            if auto %accept empty input to atomatically detect delay
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
                Im1=imtranslate(obj.rfm_b(:,:,1),[0 -obj.QS1]);
                Im2=imtranslate(obj.rfm_b(:,:,2),[0 -obj.QS2]);
                x_corr.IW=64;
                x_corr.SW=32;
                x_corr.SZ=1;
                [Im_xcorr_sl] = XCorr2D(Im1, Im2, x_corr);
                [xc_displacement,~] = obj.FindShift(Im_xcorr_sl);
                
                dy = obj.RF.speed_of_sound/obj.acq.fs;
                T = obj.acq.ftime*1E-6; %delay between pulses (seconds)
                theta = obj.flw.theta;
                dv = (dy/(T*cos(theta)))*1E3; %mm/s                
                QS_offset=median(median(xc_displacement(1:25,1:30)))/dv;

                obj.QS2 = obj.QS2-QS_offset;
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
            obj.Y = (0:dy:(size(obj.rfm,1)-1)*dy)*1E3;
            
            close(h); 
         
            
        end 
        function Detrend(obj) 
            %Remove mean of each row (normalise rows)
            obj.rfm=bsxfun(@minus,obj.rfm,mean(obj.rfm,1));       
        end
        function Init(obj)
            obj.LoadRFM;
            obj.QSCorrect('QS1',580,'QS2',580);
            obj.KWaveInit('Upsample',2);
            obj.Detrend();
%             obj.FT(0);
%             obj.BF(0);
%             obj.Highpass(5,1);
%             obj.Wallfilter();
%             obj.PlotRFM('Filter',1)            
%             obj.EnsembleCorrelation('AllCorrelations',true);
%             obj.PlotXC('SaveFig',true);
%             obj.Save();
        end
        %Reconstruction
        function TR(obj,N) %Reconstruction via time reversal
            if ~N
                N=size(obj.rfm,3);
            end
            obj.p0_recon_TR=zeros(size(obj.rfm,1),size(obj.rfm,2),N);
            h = waitbar(0, 'Initialising Waitbar');
            msg='Computing time reversal...';
            obj.KWaveInit;
            
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
            obj.p0_recon = obj.p0_recon_TR;
            disp('Done.');
        end
        function FT(obj,N,varargin) %Reconstruction via fourier transform
            if ~N
                N=size(obj.rfm,3);
            end
            padX=1.5*size(obj.rfm,2); %zero padding outside of frame to prevent wrapping         
            padY=0;

            save_opt = 0;
            if ~isempty(varargin)
                for input_index = 1:2:length(varargin)
                    switch varargin{input_index}
                        case 'Save'
                            save_opt = varargin{input_index + 1};
                        case 'PadY'                            
                            padY = varargin{input_index + 1};
                        case 'PadX'
                            padX = varargin{input_index + 1};
                        otherwise
                            error('Unknown optional input');
                    end
                end
            end
            dx = obj.kgrid.dx; 
            obj.p0_recon_FT=zeros(size(obj.rfm,1),size(obj.rfm,2),N);
            h = waitbar(0, 'Initialising Waitbar');
            msg='Computing FFT...';
            dimX=size(obj.rfm(:,:,1),2);
            dimY=size(obj.rfm(:,:,1),1);            
            sensor_data=zeros(padY+dimY,2*padX+dimX);
            
            obj.KWaveInit(size(sensor_data,2),size(sensor_data,1));
                        
            for i=1:N
                waitbar(i/N,h,msg);
                disp(i/N)    
                sensor_data(1:dimY,padX+1:padX+dimX)=obj.rfm(:,:,i);
                [recon] = kspaceLineRecon(sensor_data, dx, obj.dt, ...
                    obj.medium.sound_speed, 'Interp', '*linear');  
                
                obj.p0_recon_FT(:,:,i) = recon(1:dimY,padX+1:padX+dimX);
            end
            close(h); 
            
            if save_opt
                disp('Saving..');
                obj.Save()
                disp('Done.');
            end
            obj.p0_recon = obj.p0_recon_FT;
        end
        function BF(obj,N,varargin) %Reconstruction via beam forming
            %Code From Pim van den Berg (University of Twente)
            if ~N
                N=size(obj.rfm,3);
            end
            if ~isempty(varargin)
                aptr=varargin{1};
            else
                aptr=size(obj.rfm,2);
            end
            aptr=ceil((aptr-1)/2.0);
            
            obj.p0_recon_BF=zeros(size(obj.rfm,1),size(obj.rfm,2),N);
            
            rows=size(obj.rfm,1);
            columns=size(obj.rfm,2);                       
            pitch = obj.RF.pitch*1E-3.* obj.acq.fs/obj.RF.speed_of_sound; 
            pitch = pitch/(size(obj.rfm,2)/size(obj.rfm_b,2));

            h = waitbar(0, 'Initialising Waitbar');
            msg='Computing Beam Forming...';
            for i=1:N
                data=squeeze(obj.rfm(:,:,i));
                waitbar(i/N,h,msg);
                for channel_i = 1:columns % for every channel
                    result_xx = zeros(1,rows);
                    for dy = 1:rows % for all rows
                        for dx = 1:columns
                            if(dx ~= channel_i) % for all off_axis image points
                                if abs(channel_i - dx) <= aptr
                                delta_x = (channel_i - dx) * pitch; % Pythagoras
                                delta_y = dy;
                                % distance from imagepoint to element:
                                l = sqrt(delta_x^2 + delta_y^2);
                                l = round(l); % rounding distance
                                % to whole sample point,
                                % no interpolation
                                % check that distance is within recon area:
                                if(l >= 1 && l <= rows)
                                    result_xx(dy) = result_xx(dy) + data(l,dx); % sum image points
                                end
                                end
                            else % for the on-axis image points
                                result_xx(dy) = result_xx(dy) + data(dy, channel_i); % sum image points
                            end
                        end
                    end
                    obj.p0_recon_BF(:, channel_i, i) = result_xx;
                end
            end
            close(h);
            obj.p0_recon = obj.p0_recon_BF;
            
        end
        %Filter
        function LoadRecon(obj,varargin)
            if ~isempty(varargin)
                switch varargin{1}
                    case 'FT'
                        obj.p0_recon = obj.p0_recon_FT;
                    case 'TR'
                        obj.p0_recon = obj.p0_recon_TR;
                    case 'BF'
                        obj.p0_recon = obj.p0_recon_BF;
                    otherwise
                        error('Unknown optional input');
                end
            else
                obj.p0_recon = obj.p0_recon_FT;
            end
            
            
        end
        function Wallfilter(obj)
            obj.p0_recon(:,:,1:2:end-1)=bsxfun(@minus,...
                obj.p0_recon(:,:,1:2:end-1),mean(obj.p0_recon(:,:,1:2:end-1),3));
            obj.p0_recon(:,:,2:2:end)=bsxfun(@minus,...
                obj.p0_recon(:,:,2:2:end),mean(obj.p0_recon(:,:,2:2:end),3));
        end
        function Highpass(obj,Fp,Fst) %enter in MHz: Passband, Stopband
            %pass band frequency in rad/sample
            Fp = Fp*1E6;
            Fst = Fst*1E6;            
            hd = designfilt('highpassiir','StopbandFrequency',Fst, ...
                'PassbandFrequency',Fp,'PassbandRipple',0.2, ...
                'SampleRate',obj.acq.fs);
            %Apply filter
            obj.p0_recon = filtfilt(hd,obj.p0_recon);
        end
        %Correlation
        function EnsembleCorrelation(obj,varargin)
            all_corrs = 0;
            if isempty(varargin)
                im_stack = obj.p0_recon;
            else
                for input_index = 1:2:length(varargin)
                switch varargin{1}
                    case 'FT'
                        im_stack=obj.p0_recon_FT;
                    case 'TR'
                        im_stack=obj.p0_recon_TR;
                    case 'Raw'
                        im_stack=obj.rfm;
                    case 'AllCorrelations'
                        if ~exist('im_stack','var')
                            im_stack = obj.p0_recon;
                        end
                        all_corrs = varargin{input_index + 1};
                    otherwise
                        error ('Unkown Type: Select FT or TR')
                end
                end
            end
            assert(rem(size(im_stack,3),2)==0);
            
            x_corr.IW=32;
            x_corr.SW=16;
            x_corr.SZ=1;
            
            temp=XCorr2D(im_stack(:,:,1),im_stack(:,:,2),x_corr);
            n_corrs=size(im_stack,3)/2;
            xc_stack=zeros([size(temp),n_corrs]);
            
            %Create Time basis for xcorr
            dy = x_corr.SZ*obj.RF.speed_of_sound/obj.acq.fs;
            Y0 = x_corr.IW/2*dy;
            Yend = Y0+(size(temp,3)-1)*dy;
            obj.Y_xc = (Y0:dy:Yend)*1E3;
            
            h = waitbar(0, 'Initialising Waitbar');
            msg='Calculating Cross-Correlations...';
            for i = 1:2:size(im_stack,3)-1
                waitbar(i/(size(im_stack,3)-1),h,msg);
                xc = XCorr2D(im_stack(:,:,i),im_stack(:,:,i+1),x_corr);
                xc_stack(:,:,:,(i+1)/2) = xc;
                if all_corrs
                    [obj.xc_disp, obj.xc_amp] = obj.FindShift(xc);
                    %deconvolve xcorr amplitude
                    kernel=ones(x_corr.IW,1);
                    obj.xc_amp = deconvlucy(obj.xc_amp,kernel);
                    obj.FindFlow();
                    fig_name = ['PA_xcorr',num2str((i+1)/2)];
                    obj.PlotXC('SaveFig',true, 'FigName', fig_name);
                end
            end
            close(h);    
            %ensemble correlations:
            obj.xc_raw=squeeze(mean(xc_stack,4));
            disp(size(xc_stack));
            M = max(xc_stack,[],1);
            M = squeeze(M);
            me = median(M,3);
            M = min(M,[],3);
            
            obj.FindShift();
            obj.xc_amp = M'; %minimum intensity... 
            obj.xc_amp_me = me'; %or median projection
            
            %deconvolve xcorr amplitude
            kernel=ones(x_corr.IW,1);
            obj.xc_amp = deconvlucy(obj.xc_amp,kernel);            
            
            %estimate flow speed
            obj.FindFlow();
        end   
        %Graphical output
        function PlotRFM (obj, varargin)
            
            SaveFig = 0;
            Average = 0;
            figname = 0;
            filt = 0;
            
            if ~isempty(varargin)
                for input_index = 1:2:length(varargin)
                    switch varargin{input_index}
                        case 'SaveFig'
                            SaveFig = varargin{input_index + 1};
                        case 'FigName'
                            figname = varargin{input_index + 1};
                        case 'Average'
                            Average = varargin{input_index + 1};
                        case 'Filter'
                            filt = varargin{input_index + 1};
                        otherwise
                            error('Unknown optional input');
                    end
                end
            end
            
            if Average
                Im1=mean(obj.rfm,3);
                Im2=mean(obj.p0_recon_TR,3);
                title1='Avg Time Reversal';
                Im3=mean(obj.p0_recon_FT,3);
                title2='Avg Fourier Transform';
            else                
                Im1=obj.p0_recon_TR(:,:,1);   
                if filt
                    Im2=obj.p0_recon_FT(:,:,1);
                    title1='Unfiltered Reconstruction';
                    Im3=obj.p0_recon(:,:,1);
                    title2='Filtered Reconstruction';
                else
                    Im2=obj.p0_recon_BF(:,:,1);
                    title1='Beam Forming';
                    Im3=obj.p0_recon_FT(:,:,1);
                    title2='Fourier Transform';
                end
            end
            
            fig=figure('Visible','off');
            colormap('gray');
            
            subplot(1,3,1)
            imagesc(obj.X,obj.Y,Im1);
            title('Time Reversal');
            caxis([-150 150])
            xlabel('Lateral (mm)');
            ylabel('Depth (mm)');            
            
            subplot(1,3,2)
            imagesc(obj.X,obj.Y,Im2);
            title(title1);
            caxis([-40 40])
            caxis([-3000 3000])
            xlabel('Lateral (mm)');
            ylabel('Depth (mm)');            
            
            subplot(1,3,3)
            imagesc(obj.X,obj.Y,Im3);
            title(title2);            
            if ~filt
            caxis([-80 80])
            else
            caxis([-40 40])
            end

            xlabel('Lateral (mm)');
            ylabel('Depth (mm)');
            
            set(fig, 'Position', [100 100 800 600]);
            
            if SaveFig
                if ~figname
                    figname = [obj.pathname,'/',obj.filename,'RFM.png'];
                else
                    figname = [obj.pathname,'/',figname,'.png'];
                end
                set(gcf,'PaperPositionMode','auto')
                print(fig,figname,'-dpng','-r0')
                close(fig);
            else
                fig.Visible='on';
            end
        end
        function PlotXC (obj, varargin)
            
            SaveFig = 0;
            figname = 0;
            mask = 1;
            if ~isempty(varargin)
                for input_index = 1:2:length(varargin)
                    switch varargin{input_index}
                        case 'SaveFig'
                            SaveFig = varargin{input_index + 1};
                        case 'FigName'
                            figname = varargin{input_index + 1};
                        case 'Mask'
                            mask = varargin{input_index + 1};
                        otherwise
                            error('Unknown optional input');
                    end
                end
            end
            
            Im2=obj.xc_disp;
            Im3=obj.xc_amp;
            
            %find flowrate in mm/s
            d=obj.flw.tube_diameter; %(mm)
            A=pi*(d/2)^2;
            v=obj.flw_r/A*1000/60;
            
            fig=figure('Visible','off');            
            
            sb1 = subplot(1,2+mask,1);
            Im2=-Im2; %minus sign because tube was directed the wrong way
            imagesc(obj.X,obj.Y_xc,Im2);
%             title('Flow Speed (normalised)')
            title(['Flow Speed: ',num2str(v),'mm/s']);
            xlabel('Lateral (mm)');
            ylabel('Depth (mm)');
            load('cm_surf.mat');
            caxis([-4*v 4*v])
%             caxis([-2 2]);
            colormap(sb1,cm_surf);
            c=colorbar;
            ylabel(c,'Flow Speed (mm/s)');

            ylim([4 9]);xlim([-2 2]);
            
            sb2 = subplot(1,2+mask,2);
            imagesc(obj.X,obj.Y_xc,Im3);
            title('X-Corr Amplitude (a.u.)');
            xlabel('Lateral (mm)');
            ylabel('Depth (mm)');
            colormap(sb2,'gray');
            cmin=min(min(Im3(50:end,:)));
            cmax=max(max(Im3(50:end,:)));
            caxis([cmin cmax])
            colorbar;
            ylim([4 9]);xlim([-2 2]);
            
            if mask
                Im4 = -1*obj.xc_mask.*obj.xc_disp;
                sb3 = subplot(1,2+mask,3);
                imagesc(obj.X,obj.Y_xc,Im4)
                title(['Masked estimate: ',num2str(obj.xc_flw),'mm/s']);
                xlabel('Lateral (mm)');
                ylabel('Depth (mm)');
                caxis([-4*v 4*v])
                colormap(sb3,cm_surf);
                c=colorbar;
                ylim([4 9]);xlim([-2 2]);
                ylabel(c,'Flow Speed (mm/s)');
            end
            
            
            set(fig, 'Position', [100 100 600+400*mask 400]);
            
            if SaveFig
                if ~figname
                    figname = [obj.pathname,'/',obj.filename,'_xcorr.png'];
                else
                    figname = [obj.pathname,'/',figname,'.png'];
                end
                set(gcf,'PaperPositionMode','auto')
                print(fig,figname,'-dpng','-r0')
                close(fig);
            else
                fig.Visible='on';
            end
        end
        function PlotRecon (obj,i, varargin)
            SaveFig = 0;
            figname = 0;
            if ~isempty(varargin)
                for input_index = 1:2:length(varargin)
                    switch varargin{input_index}
                        case 'SaveFig'
                            SaveFig = varargin{input_index + 1};
                        case 'FigName'
                            figname = varargin{input_index + 1};
                        otherwise
                            error('Unknown optional input');
                    end
                end
            end
            
            Im1=obj.p0_recon(:,:,2*i-1);
            Im2=obj.p0_recon(:,:,2*i);
            
            
            fig1 = figure('Visible','off'); 
            subplot (1, 2, 1)
            colormap('gray');
            imagesc(obj.X,obj.Y,Im1);
            title('Image 1');
            caxis([-120 120])
            xlabel('Lateral (mm)');
            ylabel('Depth (mm)');            
            ylim([4 9]);xlim([-2 2]);
            
            subplot (1, 2, 2)
            colormap('gray');
            imagesc(obj.X,obj.Y,Im2);
            title('Image 2');
            caxis([-120 120])
            xlabel('Lateral (mm)');
            ylabel('Depth (mm)');
            ylim([4 9]);xlim([-2 2]);
            
            set(fig1, 'Position', [100 100 500 400]);
            if SaveFig
                if ~figname
                    figname = [obj.pathname,'/',obj.filename,num2str(i),'.png'];
                else
                    figname = [obj.pathname,'/',figname,'.png'];
                end
                set(gcf,'PaperPositionMode','auto')
                print(fig1,figname,'-dpng','-r0')
                close(fig1);
            else
                fig1.Visible='on';
            end
        end
        %helper methods     
        function Upsample(obj,N) %Upsample number of detectors
            if ~N
                error('Upsample multiplier needs to be a postive integer');
            end
            rfm_t = obj.rfm;
            for i=1:size(rfm_t,2)
                for j=1:N
                    %nearest neighbour interpolation
                    obj.rfm(:,N*(i-1)+j,:)=rfm_t(:,i,:);
                end
            end
        end         
        function Save(obj) %save the object to original folder
            save([obj.pathname,'/',obj.filename,'.mat'],'obj')
        end
        function [xc_displacement, xc_amplitude] = FindShift(obj, varargin)
            %Takes output from Xcorr2D and finds position and amplitude of
            %maxima
            if isempty(varargin)
                xc = obj.xc_raw;
            else
                xc = varargin{1};
            end
            %xc(xcorr>rows,columns,step)
            res = 100; %interpolate to 2 significant digits
            assert(rem(size(xc,1),2)~=0) %odd number of elements in colums
                                   
            dy = obj.RF.speed_of_sound/obj.acq.fs;
            T = obj.acq.ftime*1E-6; %delay between pulses (seconds)
            theta = obj.flw.theta;
            dv = (dy/(T*cos(theta*pi/180)))*1E3; %mm/s
            
            L = (size(xc,1)-1)/2;
            y = (-L:L)*dv;
            yi = (-L:1/res:L)*dv;
            
            disp('Interpolating Cross-Correlations...');
            xc_i = interp1(y,xc,yi,'spline');
            [M, I] = max(xc_i,[],1);      
            M = squeeze(M);
            I = squeeze(I);
            disp('Done.');
            
            if isempty(varargin)
                obj.xc_disp = yi(I)';
                obj.xc_amp = M';
            else
                xc_displacement = yi(I)';
                xc_amplitude = M';
            end


        end
        function FindFlow(obj, varargin)
            %Parse inputs
            cut = 30; %throw away first 30 pixels
            if ~isempty(varargin)
                for input_index = 1:2:length(varargin)
                    switch varargin{input_index}
                        case 'LaserNoise'
                            cut = varargin{input_index + 1}; %in mm
                            dy = obj.RF.speed_of_sound*obj.dt*1E3;
                            cut = cut/dy; %mm to pixels
                        otherwise
                            error('Unknown optional input');
                    end
                end
            end
            %find mask
            m = max(max(obj.xc_amp(cut:end,:)));
            m = m/2;
            obj.xc_mask = obj.xc_amp>m;
            obj.xc_mask(1:cut,:) = false;
            
            x = obj.xc_disp(obj.xc_mask);
            obj.xc_all = x;
            obj.xc_flw = median(x);
            obj.xc_flw_std = std(x);
            
        end
    end
end

