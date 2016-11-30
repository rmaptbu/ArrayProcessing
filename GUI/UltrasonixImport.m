function rfObj = UltrasonixImport(varargin)
if ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'Path'
                pathn = varargin{input_index + 1};
            otherwise
                error('Unknown optional input');
        end
    end
else
    pathn=uigetdir();
end
path=[pathn,'/'];

sample_frequency = 40e6; %[Hz]
sound_speed = 1490;      %[m/s]
pitch=0.100e-3;             %[m]
%  T=22;
%  medium.sound_speed=1.40238744e3+5.03836171*T-5.81172916e-2*(T.^2)+3.34638117e-4*(T^3)-...
%      1.48259672e-6*(T^4)+3.16585020e-9*(T^5);

% set the input arguements: force the PML to be outside the computational
% grid; switch off p0 smoothing within kspaceFirstOrder2D

[header, data] = readDAQ(path, ones(1,128), 2, true);
frameNo = header(2);
data1 = zeros(header(3)+160,(header(1)+1));

rfObj = RFObj(sample_frequency,sound_speed,pitch,...
    pathn,'RFobj',header(1)+1,header(3)+160);

for i=1:frameNo;
    [header, data] = readDAQ(path, ones(1,128), i, true);
    data1(161:end,:) =  data;
    data1(1:200,:) = 0;
    rfObj.popOn(data1);
end




