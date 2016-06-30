function [v, e]= IterateFolders()
pathname = uigetdir('/Users/Thore/Documents/Transducer Measurements/160707LinArray/');
cases=struct2cell(dir(pathname));
cases=cases(1,3:end);
v=[];
e=[];
for i=1:length(cases)
    path=['pathname/',cases{i}];
    if exist(path,'dir')
        disp(cases{i});
        [PAFrames, ~, ~] = LoadFiles(path);
        PAFrames.Init
        v = [v,PAFrames.xc_flow];
        e = [e,PAFrames.xc_flow_std];
    end
end
end