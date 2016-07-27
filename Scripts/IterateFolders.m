function [v, e]= IterateFolders()
pathname = uigetdir('/Users/Thore/Documents/Transducer Measurements/160707LinArray/');
cases=struct2cell(dir(pathname));
cases=cases(1,3:end);
v=[];
e=[];
for i=1:length(cases)
    path=[pathname,'/',cases{i}];
    if exist(path,'dir')
        disp(cases{i});
        [PAFrames, ~, ~] = LoadFiles('Path',path,'LoadExisting',1);
        if isempty(PAFrames.xc_flw)
            PAFrames.Init
        end
        v = [v,PAFrames.xc_flw];
        e = [e,PAFrames.xc_flw_std];
    end
end
end