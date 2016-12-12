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
        function svd(obj)
            for i = 1:size(obj.rfm,3);                
                sensor_data_noise=obj.rfm(:,:,i)';
                PCA=1;
                %From Xia Wenfeng (UCL)
                mean_vals=mean(sensor_data_noise, 2);
                mean_vals_mat=repmat(mean_vals',size(sensor_data_noise,2),1);
                sensor_data_noise=sensor_data_noise-mean_vals_mat';
                %
                [U,S,V] = svds(sensor_data_noise', 50);
                U=U(:,1:PCA);
                invMat=inv(U'*U)*U';
                sensor_data_t=sensor_data_noise';
                
                n_elements=size(sensor_data_t,2);
                n_time=size(sensor_data_t,1);
                
                % new filtered image
                sensor_data_filtered=zeros(n_time,n_elements);
                
                for k=1:n_elements
                    c=invMat*sensor_data_t(:,k);
                    sensor_data_filtered(:,k)=sensor_data_t(:,k)-U*c;
                end
                sensor_data=sensor_data_filtered';
                sensor_data_PCA(:,:,PCA) = sensor_data;
                obj.rfm(:,:,i)=sensor_data_PCA';
            end            
        end
        function save(obj) %save the object to original folder
            save([obj.pathname,'/',obj.filename,'.mat'],'obj')
        end
        function plot(obj,varargin)
            if ~isempty(varargin)
                fig = varargin{1};
                ax = varargin{2};
                fig.CurrentAxes=ax;
            else
                fig = figure;
            end
            Im1=obj.rfm(:,:,2);
            Im1=abs(hilbert(Im1)); %get envelope
            Im1=20*log(Im1/max(Im1(:))); %dB scale
            colormap('gray');
            imagesc(obj.X,obj.Y,Im1);
            title('Raw');
            caxis([-80 -30])
            ylim([obj.del inf])
            xlabel('Lateral (mm)');
            ylabel('Depth (mm)');            
        end
    end
end

