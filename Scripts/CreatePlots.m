pathname = uigetdir('/Users/Thore/Documents/LabData/Transducer Measurements/160707LinArray/');
cases=struct2cell(dir(pathname));
cases=cases(1,3:end);
vFT=[];
eFT=[];
vTR=[];
eTR=[];
v=[];
for i=1:length(cases)
    path=[pathname,'/',cases{i}];
    disp(i/length(cases))
    if exist(path,'dir')
        disp(cases{i});
        [PAFrames, ~, ~] = LoadFiles('Path',path,'LoadExisting',1);
        mkdir([path,'/imagesFT_filt'])
        mkdir([path,'/imagesTR_filt'])
        PAFrames.pathname=path;
        d=PAFrames.flw.tube_diameter; %(mm)
        A=pi*(d/2)^2;
        v=[v,PAFrames.flw_r/A*1000/60];
        
        %FFT
        PAFrames.LoadRecon('FT');
        PAFrames.Highpass(5,1);
        PAFrames.Wallfilter;
        PAFrames.EnsembleCorrelation;
        PAFrames.PlotXC('SaveFig',true, 'FigName', 'PAF_xcorr_FFT_filt_min');
        PAFrames.pathname=[path,'/imagesFT_filt']; 
%         for i=1:size(PAFrames.rfm,3)/2
%             PAFrames.PlotRecon(i,'SaveFig',true);
%         end
        vFT = [vFT,PAFrames.xc_flw];
        eFT = [eFT,PAFrames.xc_flw_std];
        
        %TR
        PAFrames.pathname=path;
        PAFrames.LoadRecon('TR');
        PAFrames.Highpass(5,1);
        PAFrames.Wallfilter;
        PAFrames.EnsembleCorrelation;
        PAFrames.PlotXC('SaveFig',true, 'FigName', 'PAF_xcorr_TR_filt_min');
        PAFrames.pathname=[path,'/imagesTR_filt'];
%         for i=1:size(PAFrames.rfm,3)/2
%             PAFrames.PlotRecon(i,'SaveFig',true);
%         end
        vTR = [vTR,PAFrames.xc_flw];
        eTR = [eTR,PAFrames.xc_flw_std];
        
        PAFrames.pathname=path;
    end
end
save([pathname,'/meas_filt_min.mat'],'v','vFT','eFT','vTR','eTR')
