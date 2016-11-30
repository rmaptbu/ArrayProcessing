classdef XCorrObj < handle
      properties
        %settings
        fs %sampling frequency
        c %speed of sound
        pitch %pitch between detectors
        pathname %location of data on disk
        filename %name of object when saved
        
        %data
        xc_disp %cross correlation displacement map
        xc_amp %cross correlation amplitude map
        xc_mask
        flw

        %Plots
        X %Ticks for X-axis scaling
        Y %Ticks for Y-axis scaling
        crop = [-inf inf]
    end

    methods
        function obj = XCorrObj (fs, c, pitch, pathname, filename, X, Y,...
                xc_disp, xc_amp, xc_mask,flw)
            obj.fs = fs;%sampling frequency (Hz)
            obj.c = c;%speed of sound (m/s)
            obj.pitch = pitch;%pitch between detectors (m)
            obj.pathname = pathname; %location of data
            obj.filename = filename; %name of .mat file
            obj.X = X;
            obj.Y = Y;
            obj.xc_disp = xc_disp;
            obj.xc_amp = xc_amp;
            obj.xc_mask = xc_mask;
            obj.flw = flw;
            
        end
        %Prepare Object
        function [v]=findRate(obj)
            d=obj.flw.tube_diameter; %(mm)
            A=pi*(d/2)^2;
            v=obj.flw.rate/A*1000/60; %(mm/s)
        end
        function plot (obj, varargin)
            
            SaveFig = 0;
            figname = 0;
            mask = 1;
            if ~isempty(varargin)
                for input_index = 1:2:length(varargin)
                    if ishandle(varargin{1})
                        fig = varargin{1};
                        ax = varargin{2};
                        fig.CurrentAxes=ax;
                    else                        
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
                        fig=figure('Visible','off');
                        set(fig, 'Position', [100 100 600+400*mask 400]);
                    end                    
                end
            else
                fig=figure('Visible','off');
                set(fig, 'Position', [100 100 600+400*mask 400]);
            end
            
            Im2=obj.xc_disp;
            Im3=obj.xc_amp;
            
            %find flowrate in mm/s
            d=obj.flw.tube_diameter; %(mm)
            A=pi*(d/2)^2;
            v=obj.flw.rate/A*1000/60; 
            
            if strcmp(fig.Visible,'off') 
            sb1 = subplot(1,2+mask,1);
            Im2=-Im2; %minus sign because tube was directed the wrong way
            imagesc(obj.X,obj.Y,Im2);
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
            imagesc(obj.X,obj.Y,Im3);
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
                imagesc(obj.X,obj.Y,Im4)
                title(['Masked estimate: ',num2str(obj.flw.meas_median),'mm/s']);
                xlabel('Lateral (mm)');
                ylabel('Depth (mm)');
                caxis([-2*v 2*v])
                colormap(sb3,cm_surf);
                c=colorbar;
                ylim([obj.crop]);xlim([-2 2]);
                ylabel(c,'Flow Speed (mm/s)');
            end 
            else
            
            Im2=-Im2; %minus sign because tube was directed the wrong way
            imagesc(obj.X,obj.Y,Im2);
%             title('Flow Speed (normalised)')
            title(['Flow Speed: ',num2str(v),'mm/s']);
            xlabel('Lateral (mm)');
            ylabel('Depth (mm)');
            load('cm_surf.mat');
            caxis([-2*v 2*v])
%             caxis([-2 2]);
            colormap(ax,cm_surf);
            c=colorbar('south','Color',[0.5 0.5 0.5]);
            ylabel(c,'Flow Speed (mm/s)');

            ylim(obj.crop);xlim([-2 2]);
            
            end
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
        end
end

