pathname = uigetdir('C:\Users\LABPC_TB\Documents\TransducerMeasurements\201605\16-05-12c -- Final alignment 1064\');
cases=struct2cell(dir(pathname));
cases=cases(1,3:end);

v_meas=cell(2,1);
v_all=cell(2,1);
e_meas=cell(2,1);

v=[];

mode='aptr';

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
        if isempty(PAFrames.p0_recon_BF) 
        PAFrames.BF(0); 
        end
        if isempty(PAFrames.p0_recon_FT)
        PAFrames.FT(0); 
        end
        
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
               
        v_all{2} = [v_all{2},{PAFrames.xc_all}];
        v_meas{2} = [v_meas{2},PAFrames.xc_flw];
        e_meas{2} = [e_meas{2},PAFrames.xc_flw_std];
        

        PAFrames.Save;        
   end
end
save([pathname,'/meas',mode,'.mat'],'v','v_meas','e_meas', 'v_all')


fig=figure; hold on; box on
plot(v,v);
for i=1:length(v_all)
errorbar(v, -v_meas{i}, e_meas{i});
end
legend(leg{:});
xlabel('Set Flow Speed (mm/s)');
ylabel('Measured Flow Speed (mm/s)');

figname = [pathname,'/measPlot',mode];
set(gcf,'PaperPositionMode','auto')
print(fig,[figname,'.png'],'-dpng','-r0')
savefig(fig,[figname,'.fig'])
% close(fig);


for s=1:length(v_all)
v_map=[];
for i=1:length(v_all{s});
    for j=1:length(v_all{s}{i})
    v_map=[v_map,v(i)];
    end
end
a=cell2mat(v_all{s}');
fig=figure; hold on; box on; grid on
xlabel('Set Flow Speed (mm/s)');
ylabel('Measured Flow Speed (mm/s)');
title(leg(s));
boxplot(-a, v_map)
plot(v,'.-')
for i=1:length(v)
    text(i, v(i)+10, num2str(length(v_all{s}{i})), 'HorizontalAlignment', 'center');
end
figname = [pathname,'/BoxPlot',num2str(s),mode];
set(gcf,'PaperPositionMode','auto')
print(fig,[figname,'.png'],'-dpng','-r0')
savefig(fig,[figname,'.fig'])
close(fig);

end

