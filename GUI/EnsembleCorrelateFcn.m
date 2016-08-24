function XCorrObj = EnsembleCorrelateFcn(ReconObj,varargin) 
        %Correlation
        function EnsembleCorrelation(obj,varargin)
            all_corrs = 0;
            save_xc = false;
            file = 'Ensemble';
            if isempty(varargin)
                im_stack = obj.p0_recon;
            else
                for input_index = 1:2:length(varargin)
                switch varargin{1}
                    case 'FT'
                        im_stack=obj.p0_recon_FT;
                    case 'TR'
                        im_stack=obj.p0_recon_TR;
                    case 'BF'
                        im_stack=obj.p0_recon_BF;
                    case 'Raw'
                        im_stack=obj.rfm;
                    case 'AllCorrelations'
                        if ~exist('im_stack','var')
                            im_stack = obj.p0_recon;
                        end
                        all_corrs = varargin{input_index + 1};
                    case 'SaveCorrs'
                        if ~exist('im_stack','var')
                            im_stack = obj.p0_recon;
                        end
                        save_xc = varargin{input_index + 1};
                    case 'FileName'
                        file = varargin{input_index + 1};
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
            
            if save_xc
                ensemble=obj.xc_raw;
                xc_stack_max = squeeze(max(xc_stack,[],1));  
                mask = obj.xc_mask;
                v_map = obj.xc_disp;
                amp_map = obj.xc_amp;
                x=interp1(1:length(obj.X),obj.X,...
                    1:length(obj.X)/size(obj.xc_mask,2):length(obj.X));
                y=interp1(1:length(obj.Y_xc),obj.Y_xc,...
                    1:length(obj.Y_xc)/size(obj.xc_mask,1):length(obj.Y_xc));
                v_pts = obj.xc_all;
                v_median = obj.xc_flw;
                v_std = obj.xc_flw_std;
                save([obj.pathname,'/',file,'.mat'],'xc_stack_max',...
                    'im_stack','ensemble','mask','v_map','x','y',...
                    'v_pts','v_pts','v_std')
            end
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
                    figname = [obj.pathname,'/',obj.filename,'Del',num2str(obj.QS1),'.png'];
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
            caxis([-2*v 2*v])
%             caxis([-2 2]);
            colormap(sb1,cm_surf);
            c=colorbar;
            ylabel(c,'Flow Speed (mm/s)');

            ylim(obj.crop);xlim([-2 2]);
            
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
            ylim([obj.crop]);xlim([-2 2]);
            
            if mask
                Im4 = -1*obj.xc_mask.*obj.xc_disp;
                sb3 = subplot(1,2+mask,3);
                imagesc(obj.X,obj.Y_xc,Im4)
                title(['Masked estimate: ',num2str(obj.xc_flw),'mm/s']);
                xlabel('Lateral (mm)');
                ylabel('Depth (mm)');
                caxis([-2*v 2*v])
                colormap(sb3,cm_surf);
                c=colorbar;
                ylim([obj.crop]);xlim([-2 2]);
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
        %helper methods     
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
            
            if isempty(obj.xc_mask_manual)
                obj.xc_mask = obj.xc_amp>m;
                obj.xc_mask(1:cut,:) = false;
            else
                obj.xc_mask=zeros(size(obj.xc_amp));
                x=interp1(1:length(obj.X),obj.X,...
                    1:length(obj.X)/size(obj.xc_mask,2):length(obj.X));
                y=interp1(1:length(obj.Y_xc),obj.Y_xc,...
                    1:length(obj.Y_xc)/size(obj.xc_mask,1):length(obj.Y_xc));
                m = obj.xc_mask_manual;
                Y_m=y>m{2}(1) & y<m{2}(2);
                X_m=x>m{1}(1) & x<m{1}(2);
                obj.xc_mask(Y_m, X_m)=true;
                obj.xc_mask=logical(obj.xc_mask);
            end
            
            x = obj.xc_disp(obj.xc_mask);
            obj.xc_all = x;
            obj.xc_flw = median(x);
            obj.xc_flw_std = std(x);
            
        end
    end


