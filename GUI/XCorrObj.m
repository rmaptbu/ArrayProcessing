classdef XCorrObj < handle
      properties
        %settings
        fs %sampling frequency
        c %speed of sound
        pitch %pitch between detectors
        pathname %location of data on disk
        filename %name of object when saved
        Nx %number of detectors
        Ny %number of samples
        
        %data
        xc_raw %Raw ensemble cross correlation
        xc_raw_min
        xc_disp %cross correlation displacement map
        xc_amp %cross correlation amplitude map
        xc_amp_me
        xc_mask
        xc_flw %calculated flow speed
        xc_flw_std %standard deviation of flow speed
        xc_all %all flow speed measurements
        xc_mask_manual={}       

        %Plots
        X %Ticks for X-axis scaling
        Y %Ticks for Y-axis scaling
    end

    methods
        function obj = XCorrObj (fs, c, pitch, pathname, filename, X, Y)
            obj.fs = fs;%sampling frequency (Hz)
            obj.c = c;%speed of sound (m/s)
            obj.pitch = pitch;%pitch between detectors (m)
            obj.pathname = pathname; %location of data
            obj.filename = filename; %name of .mat file
            obj.X = X;
            obj.Y = Y;
        end
        %Prepare Object
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
end

