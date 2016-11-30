function xcorrObj = EnsembleCorrelateFcn(reconObj,flw,varargin)
%% Pars Inputs
all_corrs = 0;

x_corr.IW=32;
x_corr.SW=16;
x_corr.SZ=1;
mask_mode = 'median';
mask_manual = {};
shift_order = false;
truncate = false;

if ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'XCorr'
                x_corr = varargin{input_index + 1};
            case 'AllCorrs'
                all_corrs = varargin{input_index + 1};
            case 'MaskMode'
                mask_mode = varargin{input_index + 1};
            case 'MaskManual'
                mask_manual = varargin{input_index + 1};
            case 'ShiftOrder'
                shift_order = varargin{input_index + 1};
            case 'Truncate'
                truncate = varargin{input_index + 1};
            otherwise
                error('Unknown optional input');
        end
    end
end

if truncate
    mask=zeros(size(reconObj.p0_recon_filt));
    m = reconObj.crop;
    Y_m=reconObj.Y>m(1) & reconObj.Y<m(2);
    Ymin=find(Y_m,1,'first'); %find index of non-zero elements
    Ymax=find(Y_m,1,'last');
    im_stack = reconObj.p0_recon_filt(Ymin:Ymax,:,:);
else
    im_stack = reconObj.p0_recon_filt;
end

if (rem(size(im_stack,3),2)~=0);
    if shift_order
        im_stack=im_stack(:,:,2:end);
    else
        im_stack=im_stack(:,:,1:end-1);
    end
end

%% Ensemble Correlation
dim1=ceil(x_corr.SW/2)*2-1;
dim2=size(im_stack,2);
dim3=ceil((size(im_stack,1)-x_corr.IW)/x_corr.SZ);
dim4=size(im_stack,3)/2;
xc_stack=zeros([dim1,dim2,dim3,dim4]);


%Create Time basis for xcorr
dY = x_corr.SZ*reconObj.c/reconObj.fs;
Y0 = x_corr.IW/2*dY;
Yend = Y0+(dim3-1)*dY;
Y = (Y0:dY:Yend)*1E3;
X = reconObj.X;
clear('dim1','dim2','dim3','dim4','Y0','Yend','dY');

h = waitbar(0, 'Initialising Waitbar');
msg='Calculating Cross-Correlations...';
for i = 1:2:size(im_stack,3)-1
    waitbar(i/(size(im_stack,3)-1),h,msg);
    xc_raw = XCorr2D(im_stack(:,:,i),im_stack(:,:,i+1),x_corr);
    xc_stack(:,:,:,(i+1)/2) = xc_raw;
    if all_corrs
        [xc_disp,xc_amp] = FindShift();
        %deconvolve xcorr amplitude
        kernel=ones(x_corr.IW,1);
        xc_amp = deconvlucy(xc_amp,kernel);
        FindFlow();
        fig_name = ['PA_xcorr',num2str((i+1)/2)];
        PlotXC('SaveFig',true, 'FigName', fig_name);
    end
end
close(h);

%mean of all xcorr is ensemble correlation
xc_raw=squeeze(mean(xc_stack,4));
[xc_disp,xc_amp] = FindShift();

%% Masking
M = max(xc_stack,[],1);
M = squeeze(M);
switch mask_mode
    case 'median'
        M = median(M,3);
        xc_amp = M';
    case 'minimum'
        M = min(M,[],3);
        xc_amp = M';
    case 'ensemble'
    otherwise
        error('Unknown optional input');
end

%deconvolve xcorr amplitude
kernel=ones(x_corr.IW,1);
xc_amp = deconvlucy(xc_amp,kernel);

%estimate flow speed
FindFlow();
%% create output object
xcorrObj = XCorrObj(reconObj.fs, reconObj.c, reconObj.pitch,...
    reconObj.pathname,'xcorr', X, Y, xc_disp,xc_amp,xc_mask,flw);

%% helper functions
    function [xc_disp,xc_amp] = FindShift()
        %Takes output from Xcorr2D and finds position and amplitude
        %of maxima
        
        %xc(xcorr>rows,columns,step)
        res = 100; %interpolate to 2 significant digits
        assert(rem(size(xc_raw,1),2)~=0) %odd number of elements in colums
        
        dy = reconObj.c/reconObj.fs;
        T = 1/flw.PRF; %delay between pulses (seconds)
        theta = flw.theta;
        dv = (dy/(T*cos(theta*pi/180)))*1E3; %mm/s
        
        L = (size(xc_raw,1)-1)/2;
        y = (-L:L)*dv;
        yi = (-L:1/res:L)*dv;
        
        disp('Interpolating Cross-Correlations...');
        xc_i = interp1(y,xc_raw,yi,'spline');
        [M, I] = max(xc_i,[],1);
        M = squeeze(M);
        I = squeeze(I);
        disp('Done.');
        
        xc_disp = yi(I)';
        xc_amp = M';
    end

    function FindFlow()
        %Parse inputs
        cut = 30; %throw away first 30 pixels
        %find mask
        m = max(max(xc_amp(cut:end,:)));
        m = m/2;
        
        if isempty(mask_manual)
            xc_mask = xc_amp>m;
            xc_mask(1:cut,:) = false;
        else
            xc_mask=zeros(size(xc_amp));
            m = mask_manual;
            Y_m=Y>m{2}(1) & Y<m{2}(2);
            X_m=X>m{1}(1) & X<m{1}(2);
            xc_mask(Y_m, X_m)=true;
            xc_mask=logical(xc_mask);
        end
        
        x = xc_disp(xc_mask);
        flw.meas_all = x;
        flw.meas_median = median(x);
        flw.meas_std = std(x);        
    end
end


