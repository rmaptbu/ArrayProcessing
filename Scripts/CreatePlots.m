pathname = uigetdir('C:\Users\LABPC_TB\Documents\TransducerMeasurements\201605\16-05-12c -- Final alignment 1064\');
cases=struct2cell(dir(pathname));
cases=cases(1,3:end);
v_meas=cell(4,1);
v_all=cell(4,1);
e_meas=cell(4,1);
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
        currentpath=pwd;cd(path);[~, flow_rate] = getFlowSettings();
        cd(currentpath);
        v=[v,flow_rate/A*1000/60];
        
        PAFrames.KWaveInit;
        PAFrames.BF(0); 
        PAFrames.FT(0); 
        PAFrames.PlotRFM('SaveFig',true);

        if isempty(PAFrames.p0_recon_TR)
        PAFrames.TR(0); 
        end
        
        PAFrames.LoadRecon('FT');
        PAFrames.Highpass(5,1);
        PAFrames.Wallfilter;
        PAFrames.EnsembleCorrelation;
        PAFrames.PlotXC('SaveFig',true,'FigName','XC_FT');       
        v_all{1} = [v_all{1},{PAFrames.xc_all}];
        v_meas{1} = [v_meas{1},PAFrames.xc_flw];
        e_meas{1} = [e_meas{1},PAFrames.xc_flw_std];  
        
        PAFrames.LoadRecon('TR'); 
        PAFrames.Highpass(5,1);
        PAFrames.Wallfilter;
        PAFrames.EnsembleCorrelation;
        PAFrames.PlotXC('SaveFig',true,'FigName','XC_TR');  
        v_all{2} = [v_all{2},{PAFrames.xc_all}];
        v_meas{2} = [v_meas{2},PAFrames.xc_flw];
        e_meas{2} = [e_meas{2},PAFrames.xc_flw_std];
        
        PAFrames.LoadRecon('BF');
        PAFrames.Highpass(5,1);
        PAFrames.Wallfilter;
        PAFrames.EnsembleCorrelation;
        PAFrames.PlotXC('SaveFig',true,'FigName','XC_BF');         
        v_all{3} = [v_all{3},{PAFrames.xc_all}];
        v_meas{3} = [v_meas{3},PAFrames.xc_flw];
        e_meas{3} = [e_meas{3},PAFrames.xc_flw_std];  
        
        PAFrames.p0_recon = PAFrames.rfm;
        PAFrames.Highpass(5,1);
        PAFrames.Wallfilter;
        PAFrames.EnsembleCorrelation;
        PAFrames.PlotXC('SaveFig',true,'FigName','XC_RAW');         
        v_all{4} = [v_all{4},{PAFrames.xc_all}];
        v_meas{4} = [v_meas{4},PAFrames.xc_flw];
        e_meas{4} = [e_meas{4},PAFrames.xc_flw_std];  
        
        PAFrames.Save;
        
    end
end
save([pathname,'/meas2.mat'],'v','v_meas','e_meas', 'v_all')

leg={'Known', 'FT', 'TR', 'BF', 'Raw'};

figure; hold on; box on
plot(v,v);
errorbar(v, -v_meas{1}, e_meas{1});
errorbar(v, -v_meas{2}, e_meas{2});
errorbar(v, -v_meas{3}, e_meas{3});
errorbar(v, -v_meas{4}, e_meas{4});
legend(leg{:});
xlabel('Set Flow Speed (mm/s)');
ylabel('Measured Flow Speed (mm/s)');

for s=1:length(v_all)
v_map=[];
for i=1:length(v_all{s});
    for j=1:length(v_all{s}{i})
    v_map=[v_map,v(i)];
    end
end
a=cell2mat(v_all{s}');
figure; hold on; box on; grid on
xlabel('Set Flow Speed (mm/s)');
ylabel('Measured Flow Speed (mm/s)');
title(leg(s+1));
boxplot(-a, v_map)
plot(-v,'.-')
for i=1:length(v)
    text(i, -v(i)+10, num2str(length(v_all{s}{i})), 'HorizontalAlignment', 'center');
end
end

