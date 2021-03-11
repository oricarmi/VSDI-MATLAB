%% main
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% generate individual maps %%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
% fname = "C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\2021.01.18\m210118.mat"; % path to .mat file (mYYMMDD)
% fname="D:\dataForComparison\181218\m181218.mat";
% fname ="H:\2021.02.17\m210217.mat";
% fname = "C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\2021.03.03 - target flankers\m210303.mat";
fname = "C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\2021.03.09 - target flankers\m210309.mat";
% fname ="C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\2021.01.18 - target flankers\m210118.mat";
addpath("C:\Users\orica\Dropbox\fcns_and_decript");
addpath('C:\Users\orica\Dropbox\master degree\codes');
addpath("C:\Users\Ori\Dropbox\fcns_and_decript")
addpath("C:\Users\Ori\Dropbox\master degree\codes")
global bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0 signal2_smooth basis params
%%
what = 9; % type of experiment. 8 - loc 8. 9 - loc 9 grid. 92 - 9 moving bars/diagonal dots, 2Hz
n=4; % experiment number 
% cd('C:\Users\Ori\Desktop\Ori\2nd degree\mtdt');
[cf1 cfn trsh0]=strt_up(fname, n); % parse .mat file (mYYMMDD) (legacy code)
Z = [];
switch what
    case 92  % moving bars
        for i=4:length(cfn)
            Z = [Z cfn{i}];
        end
    case 9 % loc 9
        for i=2:length(cfn)
            Z = [Z cfn{i}];
        end
    otherwise % loc 8 or other
        for i=3:length(cfn)
            Z = [Z cfn{i}];
        end
end 
% if input('Choose noise freqs? [0/1]')
%     ChooseNoiseFreqs(Z);
% end
defineParameters("C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\VSDI-MATLAB\paramsori.csv",what,rshp(Z)); % lab pc
ZZ = preProcess(Z); % preprocessing
implay(rshp(ZZ),20); % play video
[ZZZ,ZZZZ,ZZZZZ,beta] = GLM_VSDI(ZZ,[0.78 3.3 6.6],params.experiment.theoreticalSigs'); % glm denoising
params.experiment.ZZ = ZZ; params.experiment.ZZZ = ZZZ; % update params
mapTSCA = cell(params.experiment.N,1); % preallocate data for TSCA (after GLM denoising) maps
noise1.time = eye(params.experiment.T)/params.experiment.T; % autocorrelation matrix of white noise
% <---- perform TSCA
for i=1:params.experiment.N % iterate all conditions
    signal.time = params.experiment.theoreticalSigs(i,:); % this stimulus theoretical signal
    [projected,components,~,~,output] = tscaFunc(ZZZ,signal,noise1,params.TSCA.gammas(1:2),params.TSCA.numProj,params.TSCA.reduceComp); % call TSCA with this signal's parameters
    [~,I] = max(corr(abs(projected(1:9,:)'),signal.time')); % get the index of component with highest correlation to original signal
    mapTSCA{i} = postProcess(rshp(components(:,I)));
end
mapTSCA = cat(3,mapTSCA{:});
% --->
%% RUN UNTIL HERE
%% 
[~,retinotopicMap,maxind] = retinotopicMapFromIndividualMaps(mapTSCA,1,'Retinotopic Map'); % show retinotopic map
% <--- cluster analysis
X = cell(params.experiment.N,1); % preallocate memory
disp('select regoin of interest and double click inside when finished');
figure;RR = roipoly(rshp(retinotopicMap));
maxind = double(RR).*maxind;
for k=1:params.experiment.N % iterate the maps 
    [row,col] = find(maxind==k); % get pixels of this map (condition)
    X{k} = [row,col];
end
figure;[s,h] = silhouette(cat(1,X{:}),makeClustVector(cellfun(@(x) size(x,1),X))');
[~,~,~,~,~,DBI] = ClusterSimilarity(X);
figure;boxplot(s,'labels','TSCA and GLM','symbol','');box off; title('box plot of Silhouette values for each point');ylabel('Si');