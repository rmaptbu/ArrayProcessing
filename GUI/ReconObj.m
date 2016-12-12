classdef ReconObj < handle
    
    properties
        %settings
        fs %sampling frequency
        c %speed of sound
        pitch %pitch between detectors
        pathname %location of data on disk
        filename %name of object when saved
        
        %data
        p0_recon
        p0_recon_filt
        
        %WallfilterMode
        filter_mode = 'dual';
        %dual or single. Dual treats consecutive frames
        %separately 
        
        %Plots
        X %Ticks for X-axis scaling
        Y %Ticks for Y-axis scaling
        crop_y = [-inf inf]; %Plotting crop in Y
        crop_x = [-inf inf];
    end
    
    methods
        function obj = ReconObj (fs, c, pitch, pathname, filename, X, Y,varargin)
            obj.fs = fs;%sampling frequency (Hz)
            obj.c = c;%speed of sound (m/s)
            obj.pitch = pitch;%pitch between detectors (m)
            obj.pathname = pathname; %location of data
            obj.filename = filename; %name of .mat file
            obj.X = X;
            obj.Y = Y;
            if ~isempty(varargin)
                obj.p0_recon=varargin{1};
                obj.p0_recon_filt=varargin{1};
            end
        end
        function popOn (obj, p0_recon)
            if ~isempty(obj.p0_recon)
                obj.p0_recon  = cat(3,obj.p0_recon,p0_recon );
            else
                obj.p0_recon  = p0_recon ;
            end
            obj.p0_recon_filt = obj.p0_recon;
        end %Read RFM files
        function popOff (obj, N)
            assert(N<=size(obj.p0_recon,3))
            obj.p0_recon  = obj.p0_recon (:,:,1+N:end);
        end %remove N RFM files
        %Filter
        function recoverUnfiltered(obj)
            obj.p0_recon_filt=obj.p0_recon;
        end
        function wallfilter(obj)
            switch obj.filter_mode
                case 'single'
                    obj.p0_recon_filt(:,:,1:end)=bsxfun(@minus,...
                        obj.p0_recon_filt(:,:,1:end),mean(obj.p0_recon_filt(:,:,1:end),3));
                case 'dual'
                    obj.p0_recon_filt(:,:,1:2:end-1)=bsxfun(@minus,...
                        obj.p0_recon_filt(:,:,1:2:end-1),mean(obj.p0_recon_filt(:,:,1:2:end-1),3));
                    obj.p0_recon_filt(:,:,2:2:end)=bsxfun(@minus,...
                        obj.p0_recon_filt(:,:,2:2:end),mean(obj.p0_recon_filt(:,:,2:2:end),3));
                otherwise
                    error('obj.filter_mode not valid. Define dual or single');
            end
        end
        function highpass(obj,Fp,Fst) %enter in MHz: Passband, Stopband
            %pass band frequency in rad/sample
            %filters along fast time axis
            if isempty(obj.p0_recon_filt)
                obj.p0_recon_filt = obj.p0_recon;
            end
            Fp = Fp*1E6;
            Fst = Fst*1E6;
            hd = designfilt('highpassiir','StopbandFrequency',Fst, ...
                'PassbandFrequency',Fp,'PassbandRipple',0.2, ...
                'SampleRate',obj.fs);
            %Apply filter
            obj.p0_recon_filt = filtfilt(hd,obj.p0_recon_filt);
        end
        function lowpass_h(obj) %lowpass horizontally
            if isempty(obj.p0_recon_filt)
                obj.p0_recon_filt = obj.p0_recon;
            end
            hd = designfilt('lowpassiir','StopbandFrequency',0.3, ...
                'PassbandFrequency',0.1,'PassbandRipple',0.2, ...
                'SampleRate',1);
            %Apply filter
            %swap axis so filter is applied horizontally
            filt_data = permute(obj.p0_recon_filt, [2, 1, 3]);
            filt_data = filtfilt(hd,filt_data);
            obj.p0_recon_filt = permute(filt_data, [2, 1, 3]);
        end
        function cropRecon(obj,Xcrop, Ycrop) %Remove useless data
            ixmin=find((Xcrop(1)<obj.X),1,'first')-1;
            ixmax=find((Xcrop(2)>obj.X),1,'last')+1;
            iymin=find((Ycrop(1)<obj.Y),1,'first')-1;
            iymax=find((Ycrop(2)>obj.Y),1,'last')+1;
            obj.X=obj.X(ixmin:ixmax);
            obj.Y=obj.Y(iymin:iymax);
            obj.p0_recon=obj.p0_recon(iymin:iymax,ixmin:ixmax,:);
            obj.p0_recon_filt=obj.p0_recon_filt(iymin:iymax,ixmin:ixmax,:);
        end
        %Correlation
        %Graphical output
        function plot(obj,varargin)
            i = 1;
            SaveFig = 0;
            figname = 0;
            env = false;
            phase = false;
            %parse arguments
            if ~isempty(varargin)
                for input_index = 1:2:length(varargin)
                    if ishandle(varargin{1});
                        fig = varargin{1};
                        ax = varargin{2};
                        fig.CurrentAxes=ax;
                        break
                    else
                        switch varargin{input_index}
                            case 'SaveFig'
                                SaveFig = varargin{input_index + 1};
                            case 'FigName'
                                figname = varargin{input_index + 1};
                            case 'Index'
                                i = varargin{input_index + 1};
                            case 'Envelope'
                                env = varargin{input_index + 1};
                            case 'Phase'
                                phase = varargin{input_index +1};
                            otherwise
                                error('Unknown optional input');
                        end
                        fig = figure('Visible','off');
                    end
                end
            else
                fig = figure('Visible','off');
            end            
            colormap('gray');
            %find region of interest in pixels
            ixmin=find((obj.crop_x(1)<obj.X),1,'first');
            ixmax=find((obj.crop_x(2)>obj.X),1,'last');
            iymin=find((obj.crop_y(1)<obj.Y),1,'first');
            iymax=find((obj.crop_y(2)>obj.Y),1,'last');
            
            Im1=(obj.p0_recon_filt(:,:,i));
            if phase
                if env
                    error('Cannot display envelop and phase simultaneously');
                end
                Im1=hilbert(Im1); %get envelope (complex, columnwise)                
                phasemap;
                alpha=abs(Im1)/max(max(abs(Im1(iymin:iymax,ixmin:ixmax))));
                alpha=20*log(alpha); %dB mapping
                alpha=alpha+60;
                alpha(alpha<0)=0;
                alpha=alpha/60;
                Im1=angle(Im1);
            elseif env
                Im1=hilbert(Im1); %get envelope (complex, columnwise)
            end
            Im1=20*log(abs(Im1)/max(max(abs(Im1(iymin:iymax,ixmin:ixmax))))); %dB scale
            fig.Renderer = 'painters'; %faster rendering in 2D
                
            if phase
                black=cat(3,zeros(size(Im1)),zeros(size(Im1)),zeros(size(Im1)));
                image(obj.X,obj.Y, black);
                hold on
                h=imagesc(obj.X,obj.Y,Im1);
                set(h,'AlphaData',alpha);
            else
                h=imagesc(obj.X,obj.Y,Im1);
            end
            title('Reconstruction');
            
            caxis([-60 0])
%             caxis([-500 500])
            axis equal tight
            ylim([obj.crop_y]);xlim([obj.crop_x]);
            xlabel('Lateral (mm)');
            ylabel('Depth (mm)');
            
            
            if SaveFig
                if ~figname
                    figname = [obj.pathname,'/',obj.filename,num2str(i),'.png'];
                else
                    figname = [obj.pathname,'/',figname,'.png'];
                end
                set(gcf,'PaperPositionMode','auto')
                print(fig,figname,'-dpng','-r100')
                close(fig);
            else
                fig.Visible='on';
            end
        end
        function plot2 (obj,i, varargin)
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
            
            Im1=obj.p0_recon_filt(:,:,2*i-1);
            Im2=obj.p0_recon_filt(:,:,2*i);
            %             Im1=obj.p0_recon_filt(:,:,i);
            
            
            fig1 = figure('Visible','off');
            %             fig1=figure;
            subplot (1, 2, 1)
            colormap('gray');
            imagesc(obj.X,obj.Y,Im1);
            title('Image 1');
            caxis([-50 50])
            xlabel('Lateral (mm)');
            ylabel('Depth (mm)');
            ylim([obj.crop_y]);xlim([obj.crop_x]);
            subplot (1, 2, 2)
            colormap('gray');
            imagesc(obj.X,obj.Y,Im2);
            title('Image 2');
            caxis([-50 50])
            xlabel('Lateral (mm)');
            ylabel('Depth (mm)');
            %             ylim([obj.crop]);xlim([-2 2]);
            
            set(fig1, 'Position', [100 100 1600 800]);
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
        function save(obj) %save the object to original folder
            save([obj.pathname,'/',obj.filename,'.mat'],'obj')
        end
    end
end
