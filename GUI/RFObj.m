classdef RFObj < handle    
    properties
        %settings
        fs %sampling frequency 
        c %speed of sound
        pitch %pitch between detectors
        pathname %location of data on disk
        filename %name of object when saved
        
        %data
        rfm %rf matrix
        del %delay before recon     
       
        %derived
        X %Ticks for X-axis scaling
        Y %Ticks for Y-axis scaling

    end

    methods
        function obj = RFObj (fs, c, pitch, pathname, filename, Nx, Ny)
            obj.fs = fs;%sampling frequency (Hz)
            obj.c = c;%speed of sound (m/s) 
            obj.pitch = pitch;%pitch between detectors (m)
            obj.pathname = pathname; %location of data
            obj.filename = filename; %name of .mat file
            
            dx = pitch;
            dy = c/fs;
            obj.X = (0:dx:(Nx-1)*dx)*1E3-(Nx/2)*dx*1E3;
            obj.Y = (0:dy:(Ny-1)*dy)*1E3;
            
            obj.del=0;
        end
        %Prepare Object
        function popOn (obj, rfm)
            if ~isempty(obj.rfm)
                obj.rfm = cat(3,obj.rfm,rfm);
            else
                obj.rfm = rfm;
            end
        end %Read RFM files
        function popOff (obj, N)
            assert(N<=size(obj.rfm,3))
            obj.rfm = obj.rfm(:,:,1+N:end);
        end %remove N RFM files
        function detrend(obj)
            %Remove mean of each row (normalise rows)
            obj.rfm=bsxfun(@minus,obj.rfm,mean(obj.rfm,1));
        end
        function save(obj) %save the object to original folder
            save([obj.pathname,'/',obj.filename,'.mat'],'obj')
        end
        function plot(obj)
            Im1=obj.rfm(:,:,1);
            figure;
            colormap('gray');
            imagesc(obj.X,obj.Y,Im1);
            title('Slice 1');
            caxis([-120 120])
            xlabel('Lateral (mm)');
            ylabel('Depth (mm)');            
        end
    end
end

