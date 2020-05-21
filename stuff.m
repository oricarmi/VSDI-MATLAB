clear;
addpath('C:\Users\orica\Dropbox\fcns_and_decript');
addpath('C:\Users\orica\Dropbox\master degree\codes');
path = 'C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\comparison results';
files = dir(path);
for i=12:length(files)
    clear bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0 signal2_smooth basis params; close all; clc;
    global bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0 signal2_smooth basis params
    load(fullfile(files(i).folder,files(i).name));
    description = Summary.description; disp(description);
    fname = input('fname:'); % take from description
    n = input('n:');
    result = Summary.result;
    params = Summary.params;
    [cf1 cfn trsh0]=strt_up(fname, n);  
    fn = fieldnames(result);
    for k=1:length(fn) % iterate the methods 
        for j=1:size(result.(fn{k}).maps,3)
            result.(fn{k}).maps(:,:,j) = MinMaxNorm((result.(fn{k}).maps(:,:,j)));
        end
        [~,result.(fn{k}).retinotopicMap] = retinotopicMapFromIndividualMaps(result.(fn{k}).maps,5,fn{k});
    end
    [result.TSCA.performance,result.Tmax.performance,result.AOF.performance,result.Corr.performance,result.GLM.performance,result.Nadav.performance] = performanceRealData(result);
    isGood = input('is good?');
    if isGood
        Summary = struct('params',params,'result',result','description',description);
        [file,path] = uiputfile;
        save(fullfile(path,file),'Summary');
    end
end
%%
ZSig = zeros(m*m,T); % preallocate memory
for i = 1:T
    ZSig(:,i) = reshape( ...
        signals(1).time(i)*signals(1).space+signals(2).time(i)*signals(2).space+...
        signals(3).time(i)*signals(3).space+signals(4).time(i)*signals(4).space+...
        signals(5).time(i)*signals(5).space,[],1);
end
ZNoise = zeros(m*m,T); % preallocate memory
for i = 1:T
    ZNoise(:,i) = reshape( ...
        noises(1).time(i)*noises(1).space+...
        noises(2).time(i)*noises(2).space+...
        noises(3).time(i)*noises(3).space,[],1);
end
