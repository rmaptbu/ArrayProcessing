pathname = uigetdir('C:\Users\LABPC_TB\Documents\TransducerMeasurements\201605\16-05-12c -- Final alignment 1064\');
cases=struct2cell(dir(pathname));
cases=cases(1,3:end);
v_meas=cell(32,1);
v_all=cell(32,1);
e_meas=cell(32,1);
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
        
        PAFrames.Init
        leg={};
        for aptr=1:32
            aperture=aptr*4;
            PAFrames.BF(0,aperture);
            PAFrames.LoadRecon('BF');
            PAFrames.Highpass(5,1);
            PAFrames.Wallfilter;
            PAFrames.EnsembleCorrelation;
            PAFrames.PlotRecon(2,'SaveFig',true,'FigName',['Recon_BF',mode,num2str(aptr)]);
            PAFrames.PlotXC('SaveFig',true,'FigName',['XC_BF',mode,num2str(aptr)]);
            v_all{aptr} = [v_all{aptr},{PAFrames.xc_all}];
            v_meas{aptr} = [v_meas{aptr},PAFrames.xc_flw];
            e_meas{aptr} = [e_meas{aptr},PAFrames.xc_flw_std];
            leg{aptr} = num2str(aperture);
        end
        
    end
end
save([pathname,'/meas',mode,'.mat'],'v','v_meas','e_meas', 'v_all')


fig=figure; hold on; box on
plot(v,v);
for i=1:5:32
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

