function [Im_xcorr_sl] = XCorr2D(Im1, Im2, x_corr, varargin) %Xcorr Im1/2 along dim1

%Pass Arguments, setup variables
SW = x_corr.SW;
IW = x_corr.IW;
SZ = x_corr.SZ;
unbiased = 1; %default to scaling xcorr for unbiased

if ~isempty(varargin)
    for input_index = 1:2:length(varargin)
        switch varargin{input_index}
            case 'unbiased'
                unbiased = varargin{input_index + 1};
            otherwise
                error('Unknown optional input');
        end
    end
end
sw = ceil(SW/2);

%Cross Correlation--------------------------------
n_corrs=ceil((size(Im1,1)-IW)/SZ);%number of xcorrs
Im1_sl=zeros([IW, size(Im1,2), n_corrs]);
Im2_sl=zeros([IW, size(Im1,2), n_corrs]);
%create array of Interrogation windows
%Xcorr only along that line.
for i=1:n_corrs
    Im1_sl(:,:,i) = Im1((i-1)*SZ+1:(i-1)*SZ+IW,:);
    Im2_sl(:,:,i) = Im2((i-1)*SZ+1:(i-1)*SZ+IW,:);
end
%Fourier Transform all sequences
Im1_sl_fft=fft(Im1_sl,2^nextpow2(2*IW-1));
Im2_sl_fft=fft(Im2_sl,2^nextpow2(2*IW-1));
%Perform cross-correlation
Im_xcorr_sl = ifft(Im1_sl_fft.*conj(Im2_sl_fft));
Im_xcorr_sl = real(Im_xcorr_sl);
%Reorder and only keep search window
Im_xcorr_sl = [Im_xcorr_sl(end-sw+2:end,:,:);Im_xcorr_sl(1:sw,:,:)];


%slow xcorr for comparison. Should give the same result...
%             tic
%             Im_xcorr_sl=zeros([2*IW-1, size(Im1,2), n_corrs]);
%             for line=1:size(Im1,2)
%                 for i=1:SZ:(size(Im1,1)-IW)
%                     Im_xcorr_sl(:,line,((i-1)/SZ)+1)=xcorr(Im1(i:i+IW-1,line),Im2(i:i+IW-1,line));
%                 end
%             end
%             toc

if unbiased
    %Scale for bias in finite length xcorr
    maxlag=IW-1;
    lags = -maxlag:maxlag;
    scale = (IW-abs(lags))';
    scale = scale(IW-sw+1:IW+sw-1);
    scale = repmat(scale,[1,size(Im_xcorr_sl,2),size(Im_xcorr_sl,3)]);
    Im_xcorr_sl = Im_xcorr_sl./scale;
end

end