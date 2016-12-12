function rfObj = UltrasonixImport(varargin)
mode = 'PA'; %Ultrasound or Photacoustic
if ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'Path'
                pathn = varargin{input_index + 1};
            case 'Mode'
                mode = varargin{input_index + 1};
            otherwise
                error('Unknown optional input');
        end
    end
else
    pathn=uigetdir();
end
path=[pathn,'/'];


sound_speed = 1490;      %[m/s]
pitch=0.100e-3;             %[m]
%  T=22;
%  medium.sound_speed=1.40238744e3+5.03836171*T-5.81172916e-2*(T.^2)+3.34638117e-4*(T^3)-...
%      1.48259672e-6*(T^4)+3.16585020e-9*(T^5);

% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D
switch mode
    case 'PA'
        sample_frequency = 40e6; %[Hz]
        [header, data] = readDAQ(path, ones(1,128), 2, true);
        frameNo = header(2);
        data1 = zeros(header(3)+160,(header(1)+1));
        Nx = header(1)+1;
        Ny = header(3)+160;
        %first 160 are empty because of delay in acquisition
        rfObj = RFObj(sample_frequency,sound_speed,pitch,...
            pathn,'RFobj',Nx,Ny);        
        for i=1:frameNo;            
            [~, data] = readDAQ(path, ones(1,128), i, true);
            data1(161:end,:) =  data;
            data1(1:200,:) = 0;
            rfObj.popOn(data1);  
        end
    case 'US'
        sample_frequency = 80e6; %[Hz]
        if ~exist([path,'TexoData.rf'],'file');
            error([path,'TexoData.rf does not exist']);
        end
        pitch=pitch/2;
        [data, header] = RPread([path,'TexoData.rf']);
        frameNo = header.nframes;
        Ny = header.h;
        Nx = header.w-2;
        dx = pitch;
        dy = sound_speed/sample_frequency;
        data=data(:,3:end,:);
        X = (0:dx:(Nx-1)*dx)*1E3-(Nx/2)*dx*1E3;
        Y = (0:dy:(Ny-1)*dy)*1E3;
        rfObj = ReconObj(sample_frequency,sound_speed,pitch,...
            pathn,'ReconObjUS',X,Y);
        rfObj.p0_recon=data;
        rfObj.p0_recon_filt=data; 
    otherwise
        error('Unkown mode. Enter either PA or US');
end







