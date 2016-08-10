pathname = uigetdir('C:\Users\LABPC_TB\Documents\TransducerMeasurements\201605\');
cases=struct2cell(dir(pathname));
cases=cases(1,3:end);
v_meas_qs2=cell(5,1);
e_meas_qs2=cell(5,1);
v_meas_padY1=cell(5,1);
e_meas_padY1=cell(5,1);
v_meas=cell(2,1);
e_meas=cell(2,1);
v=[];

for i=1:length(cases)
    
    path=[pathname,'/',cases{i}];
    disp(i/length(cases))
    if exist(path,'dir')
        
        disp(cases{i});
        [PAFrames, ~, ~] = LoadFiles('Path',path,'LoadExisting',1);
        PAFrames.pathname=path;
        d=PAFrames.flw.tube_diameter; %(mm)
        A=pi*(d/2)^2;
        v=[v,PAFrames.flw_r/A*1000/60];
        for padY1=1:5
            mkdir([path,'/imagesFT_padY1',num2str(padY1*50)])
            PAFrames.LoadRFM;
            PAFrames.QSCorrect('QS1',580,'QS2',580);
            PAFrames.KWaveInit('Upsample',4);
            PAFrames.Detrend();
            PAFrames.FT(0,'PadY1',padY1*50);
            PAFrames.LoadRecon('FT');
            PAFrames.Highpass(5,1);
            PAFrames.Wallfilter;
            PAFrames.EnsembleCorrelation;
            PAFrames.pathname=[path,'/imagesFT_padY1',num2str(padY1*50)];
            PAFrames.PlotXC('SaveFig',true, 'FigName', ['PAF_xcorr_FFT', num2str(padY1)]);
            PAFrames.PlotRecon(1,'SaveFig',true,'FigName', ['PAF_recon', num2str(padY1)]);
            v_meas_padY1{padY1} = [v_meas_padY1{padY1},PAFrames.xc_flw];
            e_meas_padY1{padY1} = [e_meas_padY1{padY1},PAFrames.xc_flw_std]; 
            PAFrames.pathname=path;
        end
        for qs2=1:5
            mkdir([path,'/imagesFT_qs2',num2str((qs2-3)*10)])
            PAFrames.LoadRFM;
            PAFrames.QSCorrect('QS1',580,'QS2',580+(qs2-3)*10);
            PAFrames.KWaveInit('Upsample',4);
            PAFrames.Detrend();
            PAFrames.FT(0);
            PAFrames.LoadRecon('FT');
            PAFrames.Highpass(5,1);
            PAFrames.Wallfilter;
            PAFrames.EnsembleCorrelation;
            PAFrames.pathname=[path,'/imagesFT_qs2',num2str((qs2-3)*10)];
            PAFrames.PlotXC('SaveFig',true, 'FigName', ['PAF_xcorr_FFT', num2str(qs2)]);
            PAFrames.PlotRecon(1,'SaveFig',true,'FigName', ['PAF_recon', num2str(qs2)]);
            v_meas_qs2{qs2} = [v_meas_qs2{qs2},PAFrames.xc_flw];
            e_meas_qs2{qs2} = [e_meas_qs2{qs2},PAFrames.xc_flw_std]; 
            PAFrames.pathname=path;
        end
        mkdir([path,'/imagesFTTR'])
            PAFrames.LoadRFM;
            PAFrames.QSCorrect('QS1',580,'QS2',580);
            PAFrames.KWaveInit('Upsample',4);
            PAFrames.Detrend();
            PAFrames.FT(0);
            PAFrames.LoadRecon('FT');
            PAFrames.Highpass(5,1);
            PAFrames.Wallfilter;
            PAFrames.EnsembleCorrelation;            
            PAFrames.pathname=[path,'/imagesFTTR'];
            PAFrames.PlotXC('SaveFig',true, 'FigName', 'PAF_xcorr_FT');
            PAFrames.PlotRecon(1,'SaveFig',true, 'FigName', 'PAF_recon_FT');
            v_meas{1} = [v_meas{1},PAFrames.xc_flw];
            e_meas{1} = [e_meas{1},PAFrames.xc_flw_std]; 


           
            PAFrames.LoadRecon('TR');
            PAFrames.Highpass(5,1);
            PAFrames.Wallfilter;
            PAFrames.EnsembleCorrelation;
   
            PAFrames.PlotXC('SaveFig',true, 'FigName', 'PAF_xcorr_TR');
            PAFrames.PlotRecon(1,'SaveFig',true, 'FigName', 'PAF_recon_TR');
            v_meas{2} = [v_meas{2},PAFrames.xc_flw];
            e_meas{2} = [e_meas{2},PAFrames.xc_flw_std]; 
            PAFrames.pathname=path;
    end
end
save([pathname,'/meas.mat'],'v','v_meas_qs2','e_meas_qs2',...
    'v_meas_padY1','e_meas_padY1','v_meas','e_meas')
