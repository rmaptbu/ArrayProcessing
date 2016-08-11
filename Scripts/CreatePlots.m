pathname = uigetdir('/Users/Thore/Documents/LabData/Transducer Measurements/160707LinArray/');
cases=struct2cell(dir(pathname));
cases=cases(1,3:end);
v_meas=cell(2,1);
v_all=cell(2,1);
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
        currentpath=pwd;cd(path);[~, flow_rate] = getFlowSettings();
        cd(currentpath);
        v=[v,flow_rate/A*1000/60];
        
        PAFrames.Init; 
        
        PAFrames.LoadRecon('FT');
        PAFrames.Highpass(5,1);
        PAFrames.Wallfilter;
        PAFrames.EnsembleCorrelation;
               
        v_all{1} = [v_all{1},{PAFrames.xc_all}];
        v_meas{1} = [v_meas{1},PAFrames.xc_flw];
        e_meas{1} = [e_meas{1},PAFrames.xc_flw_std];  
        
        PAFrames.LoadRecon('BF'); 
        PAFrames.Highpass(5,1);
        PAFrames.Wallfilter;
        PAFrames.EnsembleCorrelation;
        PAFrames.Save;

        v_all{2} = [v_all{2},{PAFrames.xc_all}];
        v_meas{2} = [v_meas{2},PAFrames.xc_flw];
        e_meas{2} = [e_meas{2},PAFrames.xc_flw_std];
        
    end
end
save([pathname,'/meas2.mat'],'v','v_meas','e_meas', 'v_all')

figure; hold on
errorbar(v, v_meas{1}, e_meas{1});
errorbar(v, v_meas{2}, e_meas{1});
 

