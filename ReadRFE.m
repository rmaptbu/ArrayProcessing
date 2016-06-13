function [finfo, acq, tran, lm, rfm] = ReadRFE(filename, nr_load_fr, nr_skip_fr)
%[finfo,rfm] = Read_RFe(filename, begin_frame, nr_frames)   Read plane wave *.rfe file   
%
%   Input
%   filename    : name of the <*.rfe> file (if 0 file dialog will open)
%   nr_load_fr and nr_skip_fr: load nr_load_fr then skip nr_skip_fr and repeat.
%
%   Output
%   finfo       : File info
%   acq         : Acquisition settings
%   tran        : Transmission settings
%   lm          : Label matrix (16,frame_nr)
%   rfm         : RF-data matrix (nr_sampels, line_nr);
%
%
%   Peter Brands, Yannick van Bavel
%	Modified by Pim van den Berg to select every several frames.

%
%   Get file name if there is no file input available
%

%
%   Open *.zre file
%
fid = fopen(filename,'r');

%
%   Start interpertation of file
%
if (fid == -1)
   error(1,'File "%s" does NOT exist in this directory!\n\n',filename)
else
    
%
%   Read file header
%
    finfo.id  = fread(fid, 10, 'uint8=>char')';     % File ID "planewave"
    finfo.ver = fread(fid, 1, 'uint32');            % File version
    finfo.size= fread(fid, 1, 'uint32');            % File size in bytes
    finfo.date= fread(fid, 20, 'uint8=>char')';     % File date and time stamp
    finfo.nrf = fread(fid, 1, 'uint32');            % Number of frames
    finfo.nrl = fread(fid, 1, 'uint32');            % Number of lines
    finfo.nrs = fread(fid, 1, 'uint32');            % Number of samples
    finfo.lmp = fread(fid, 1, 'uint32');            % Label matrix offset [bytes]
    finfo.rfp = fread(fid, 1, 'uint32');            % RF-matrix offset [bytes]\
	
	nr_frames = finfo.nrf;
	begin_frame = 1;
    
%
%   Read acqusition settings
%
    acq.size = fread(fid, 1, 'uint32');             % Size of the acqusition settings struct [bytes]
    acq.ftime = fread(fid, 1, 'uint32');            % Frame time [us]
    acq.fs = fread(fid, 1, 'uint32');               % Sample frequency [Hz]
    acq.dt = fread(fid, 1, 'uint32');               % Acqusition delay [RF sampels]

%
%    Read transmission settings
%
    tran.size = fread(fid, 1, 'uint32');
    tran.halfWaveLength = fread(fid, 1, 'uint32');
    tran.pulseDuration = fread(fid, 1, 'uint32');
    tran.numberHalfWaves = fread(fid, 1, 'uint32');
    tran.phase = fread(fid, 1, 'uint32');
    tran.pulseAmplitude = fread(fid, 1, 'uint32');
    
%
%   Skip number of labels
%
    fseek(fid, 32*(begin_frame-1), 0);
    
%
%   Read label matrix
%
    lm = reshape(fread(fid, 16*nr_frames, 'int16'), 16, nr_frames);
    
%
%   Skip number of frames
%
    fs = finfo.rfp+finfo.nrl*finfo.nrs*2*(begin_frame-1);   % Pointer to RF-frame in bytes
    fseek(fid, fs, -1);

%
%   Read RF-matrix
%
    nr_bursts = floor(nr_frames / (nr_load_fr+nr_skip_fr));
    nr_post_frames = nr_bursts * nr_load_fr;
    nrrs = nr_load_fr*finfo.nrs*finfo.nrl; % number of samples per read
% 	rfm = zeros(finfo.nrs*finfo.nrl*nr_post_frames,1,'int16');
    rfm = zeros(finfo.nrs,finfo.nrl,nr_post_frames,1,'double');

	for n = 1:nr_bursts
        nri = (n-1)*nrrs; % sample index
		rfm( nri+(1:nrrs) ) = fread(fid,nrrs, 'int16');
		fseek(fid,2*finfo.nrs*finfo.nrl*(nr_skip_fr),0);
	end
%     rfm = double(reshape(rfm, finfo.nrs, finfo.nrl, nr_post_frames));
%
%   Close the file
%       
    fclose(fid);
    
end

%
%   END
%

