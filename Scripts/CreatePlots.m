pathname = uigetdir('/Users/Thore/Documents/LabData/Transducer Measurements/160707LinArray/');
cases=struct2cell(dir(pathname));
cases=cases(1,3:end);
vFT=[];
eFT=[];
vTR=[];
eTR=[];
for i=1:length(cases)
    path=[pathname,'/',cases{i}];
    disp(i/length(cases))
    if exist(path,'dir')
        disp(cases{i});
        [PAFrames, ~, ~] = LoadFiles('Path',path,'LoadExisting',1);
        mkdir([path,'/imagesFT'])
        mkdir([path,'/imagesTR'])
        PAFrames.pathname=path;
        
        %FFT
        PAFrames.p0_recon=PAFrames.p0_recon_FT;
        PAFrames.EnsembleCorrelation;
        PAFrames.PlotXC('SaveFig',true, 'FigName', 'PAF_xcorr_FFT');
        PAFrames.pathname=[path,'/imagesFT'];
        for i=1:size(PAFrames.rfm,3)/2
            PAFrames.PlotRecon(i,'SaveFig',true);
        end
        vFT = [vFT,PAFrames.xc_flw];
        eFT = [eFT,PAFrames.xc_flw_std];
        
        %TR
        PAFrames.p0_recon=PAFrames.p0_recon_TR;
        PAFrames.EnsembleCorrelation;
        PAFrames.PlotXC('SaveFig',true, 'FigName', 'PAF_xcorr_TR');
        PAFrames.pathname=[path,'/imagesTR'];
        for i=1:size(PAFrames.rfm,3)/2
            PAFrames.PlotRecon(i,'SaveFig',true);
        end
        vTR = [vTR,PAFrames.xc_flw];
        eTR = [eTR,PAFrames.xc_flw_std];
        
        PAFrames.pathname=path;
    end
end
save([path,'/meas.mat'],'vFT','eFT','vTR','eTR')
