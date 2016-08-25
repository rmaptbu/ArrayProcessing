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
        
        %Plots
        X %Ticks for X-axis scaling
        Y %Ticks for Y-axis scaling
        crop=[-inf inf];
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
        function wallfilter(obj)
            obj.p0_recon_filt(:,:,1:end)=bsxfun(@minus,...
                obj.p0_recon_filt(:,:,1:end),mean(obj.p0_recon_filt(:,:,1:end),3));
            %
            %             obj.p0_recon_filt(:,:,1:2:end-1)=bsxfun(@minus,...
            %                 obj.p0_recon_filt(:,:,1:2:end-1),mean(obj.p0_recon_filt(:,:,1:2:end-1),3));
            %             obj.p0_recon_filt(:,:,2:2:end)=bsxfun(@minus,...
            %                 obj.p0_recon_filt(:,:,2:2:end),mean(obj.p0_recon_filt(:,:,2:2:end),3));
        end
        function highpass(obj,Fp,Fst) %enter in MHz: Passband, Stopband
            %pass band frequency in rad/sample
            if isemtpy(obj.p0_recon_filt)
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
        %Correlation
        %Graphical output
        function plot (obj,i, varargin)
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
            
            
            fig1 = figure('Visible','off'); 
            subplot (1, 2, 1)
            colormap('gray');
            imagesc(obj.X,obj.Y,Im1);
            title('Image 1');
            caxis([-120 120])
            xlabel('Lateral (mm)');
            ylabel('Depth (mm)');            
            ylim([obj.crop]);xlim([-2 2]);
            
            subplot (1, 2, 2)
            colormap('gray');
            imagesc(obj.X,obj.Y,Im2);
            title('Image 2');
            caxis([-120 120])
            xlabel('Lateral (mm)');
            ylabel('Depth (mm)');
            ylim([obj.crop]);xlim([-2 2]);
            
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
        function save(obj) %save the object to original folder
            save([obj.pathname,'/',obj.filename,'.mat'],'obj')
        end
    end
end
