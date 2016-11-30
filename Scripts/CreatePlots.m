pathname = uigetdir('/Users/Thore/Documents/LabData/Transducer Measurements/160707LinArray/16-05-12b -- 1064 nm - 40p blood in PBS/001 Range of flows 65 ns low edge/data/');
cases=struct2cell(dir(pathname));
cases=cases(1,3:end);
cf=pwd;
cd(pathname);

v_meas=cell(2,1);
v_all=cell(2,1);
e_meas=cell(2,1);

v=[];

mode='';
crop=[0 6];
for i=1:length(cases)    
    path=[pathname,'/',cases{i}];
%     disp(i/length(cases))
    if exist(path,'dir')        
        disp(cases{i});
        [PAFrames, ~, ~] = LoadFiles('Path',path,'LoadExisting',1);
        PAFrames.pathname=path;%make sure pathname is correct
        
        %find set flow speed
        d=PAFrames.flw.tube_diameter; %(mm)
        A=pi*(d/2)^2;        
        currentpath=pwd;cd(path);
        [~, flow_rate] = getFlowSettings();
        cd(currentpath);        
        v=[v,flow_rate/A*1000/60];   
        
        %Setup reconstructions
        PAFrames.Init;
        if isempty(PAFrames.p0_recon_BF)
        PAFrames.BF(0); 
        end
        if isempty(PAFrames.p0_recon_FT)
        PAFrames.FT(0);         
        end
        PAFrames.Save;          
        
        %Plot Ensemble Correlations     
        PAFrames.crop=crop; %region of interest
        PAFrames.xc_mask_manual={};
        mkdir('../fig');
        PAFrames.pathname=[pathname,'/../fig'];
        
        PAFrames.PlotRFM('SaveFig',true,'FigName',['recon',num2str(i)]);
        
        PAFrames.LoadRecon('FT');
        PAFrames.Highpass(5,1);
        PAFrames.Wallfilter;
        PAFrames.pathname=path; 
        PAFrames.EnsembleCorrelation('SaveCorrs',true);        
        PAFrames.pathname=[pathname,'/../fig'];

        PAFrames.PlotXC('SaveFig',true,'FigName',['XC_FT',num2str(i)]);       
        v_all{1} = [v_all{1},{PAFrames.xc_all}];
        v_meas{1} = [v_meas{1},PAFrames.xc_flw];
        e_meas{1} = [e_meas{1},PAFrames.xc_flw_std]; 
        
                
        PAFrames.LoadRecon('BF');
        PAFrames.Highpass(5,1);
        PAFrames.Wallfilter;
        PAFrames.pathname=path; 
        PAFrames.EnsembleCorrelation('SaveCorrs',true); ;        
        PAFrames.pathname=[pathname,'/../fig'];

        PAFrames.PlotXC('SaveFig',true,'FigName',['XC_BF',num2str(i)]);  
        v_all{2} = [v_all{2},{PAFrames.xc_all}];
        v_meas{2} = [v_meas{2},PAFrames.xc_flw];
        e_meas{2} = [e_meas{2},PAFrames.xc_flw_std];
        
        %make sure pathname is correct
        PAFrames.pathname=path;      
   end
end
save([pathname,'/meas',mode,'.mat'],'v','v_meas','e_meas', 'v_all')
cd(cf);

%%Plots
leg={'Known','FT','BF'};

fig=figure; hold on; box on
[ve,i]=sort(v);
plot(ve,ve);
errorbar(ve, -v_meas{1}(i), e_meas{1}(i));
errorbar(ve, -v_meas{2}(i), e_meas{2}(i));
legend(leg{:},'Location','northwest');
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
title(leg(s+1));
boxplot(-a, v_map)
plot(v,'.-')
for i=1:length(v)
    text(i, v(i)+10, num2str(length(v_all{s}{i})), 'HorizontalAlignment', 'center');
end
figname = [pathname,'/../fig/BoxPlot',num2str(s),mode];
set(gcf,'PaperPositionMode','auto')
print(fig,[figname,'.png'],'-dpng','-r0')
savefig(fig,[figname,'.fig'])
close(fig);

end

figure;
m=floor(size(obj.rfm,2)/2);
hold on;
box on;
grid on;
title('Profile at centre line');
plot(obj.Y,obj.p0_recon_BF(:,m,1)/max(obj.p0_recon_BF(100:end,m,1)),'LineWidth',2);
% plot(obj.Y,obj.p0_recon_TR(:,m,1)/max(obj.p0_recon_TR(100:end,m,1))+1,'LineWidth',2);
plot(obj.Y,obj.p0_recon_FT(:,m,1)/max(obj.p0_recon_FT(100:end,m,1))+2,'LineWidth',2);
legend('BF','TR','FT');
xlim(crop);
xlabel('Depth (mm)');
figname = [pathname,'/../fig/profile',mode];
set(gcf,'PaperPositionMode','auto')
print(fig,[figname,'.png'],'-dpng','-r0')
savefig(fig,[figname,'.fig'])
close(fig);
