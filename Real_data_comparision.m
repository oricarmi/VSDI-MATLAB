%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERAL %%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
% fname = "G:\2020.01.21\m200121.mat"; n=4;
% fname = "E:\2018.12.18\MAT\m181218.mat"; n=2;
% fname = "E:\191119\m191119.mat"; n=2;
% fname = "E:\200121\m200121.mat"; n=4;
% fname = "E:\180904\m180904.mat"; n=4;
fname = "E:\181218\m181218.mat"; n=2;
% fname = "E:\180801\m180801.mat"; n=3;
% fname = "D:\2019.07.10\m190710.mat"; n=2;
% fname = "G:\2019.12.18\m191218.mat"; n=1;
% fname = "G:\2019.12.11\m191211.mat"; n=3;
% cd('C:\Users\Ori\Desktop\Ori\2nd degree\mtdt');
addpath("C:\Users\orica\Dropbox\fcns_and_decript");
addpath('C:\Users\orica\Dropbox\master degree\codes');
addpath("C:\Users\Ori\Dropbox\fcns_and_decript")
addpath("C:\Users\Ori\Dropbox\master degree\codes")
global bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0 signal2_smooth basis params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result = struct('TSCA',struct,'Tmax',struct,'AOF',struct,'Corr',struct,'GLM',struct,'Nadav',struct);
[result.TSCA.maps,result.Tmax.maps,result.AOF.maps,result.Corr.maps,result.GLM.maps,result.Nadav.maps,tscaBoth] = GenerateMaps(fname,n,8); % 3rd parameter: what, 8=loc8,9=loc9,92=mvngbars2hz
% <----- generate retinotopic maps from the individual maps
fn = fieldnames(result);
for i=1:length(fn) % iterate the methods 
    [~,result.(fn{i}).retinotopicMap] = retinotopicMapFromIndividualMaps(result.(fn{i}).maps,1,fn{i});
end
% ---->
[result.TSCA.performance,result.Tmax.performance,result.AOF.performance,result.Corr.performance,result.GLM.performance,result.Nadav.performance] = performanceRealData(result);
result = clusterEvaluation(result);
Summary = struct('params',params,'result',result','description','loc 4,191119 n=2'); 
%% DB index for cluster similarity
X = cell(8,1); Xraw = cell(8,1); mapp = zeros(size(brn,1),size(brn,2),8);
for i=1:8
    [row,col] = find(result.TSCA.maps(:,:,i)>prctile(reshape(result.TSCA.maps(:,:,i),[],1),99));
    mapp(:,:,i) = result.TSCA.maps(:,:,i)>prctile(reshape(result.TSCA.maps(:,:,i),[],1),99);
    X{i} = [row,col];
    Xraw{i} = result.TSCA.maps(row,col,i);
end
[R,dbs,dbsI,dbns,dbnsI,DBI] = ClusterSimilarity(X);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Show Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath("C:\Users\orica\Dropbox\fcns_and_decript");
addpath('C:\Users\orica\Dropbox\master degree\codes');
addpath("C:\Users\Ori\Dropbox\fcns_and_decript")
addpath("C:\Users\Ori\Dropbox\master degree\codes")
global bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0 signal2_smooth basis params
% description = Summary.description
% fname = "E:\191119\m191119.mat";n=2; % take from description
% fname = "E:\181218\m181218.mat"; n=3;
% fname = "E:\200121\m200121.mat"; n=5;
result = Summary2.result;
params = Summary2.params;
% [cf1 cfn trsh0]=strt_up(fname, n);  
fn = fieldnames(result);
retmapsAll = figure("name","retmapsAll");
% indmapsAll = figure("name","indmapsAll");
for i=2:length(fn)-1 % iterate the methods    
    figure(retmapsAll); subplot(2,3,i-1)
    [~,r] = retinotopicMapFromIndividualMaps(result.(fn{i}).maps,0,fn{i},92);
    imf2(r); 
    if contains(fn{i},'nadav','ignorecase',true)
        title('M.P.T');
    elseif contains(fn{i},'tsca','ignorecase',true)
        title('TSCA');
    else
        title(fn{i});
    end
%     r = plotMaps(result.(fn{i}).maps,fn{i},1);
end
% with and without GLM
figure;
subplot 121
[~,r] = retinotopicMapFromIndividualMaps(result.TSCAnoGLM.maps,0,'TSCA without GLM',92);
imf2(r); title('TSCA without GLM');
subplot 122
[~,r] = retinotopicMapFromIndividualMaps(result.TSCAwGLM.maps,0,'TSCA with GLM',92);
imf2(r);title('TSCA after GLM denoising');
