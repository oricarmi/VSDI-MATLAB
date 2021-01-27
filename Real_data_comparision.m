%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERAL %%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
% fname = "G:\2020.01.21\m200121.mat"; n=4;
% fname = "E:\2018.12.18\MAT\m181218.mat"; n=2;
% fname = "E:\191119\m191119.mat"; n=2;
% fname = "E:\200121\m200121.mat"; n=4;
% fname = "E:\180904\m180904.mat"; n=4;
% fname ="G:\comparision results 2\loc 8,181218 n=2.mat"; n=2;
% fname = "E:\181218\m181218.mat"; n=2;
% fname = "E:\180801\m180801.mat"; n=3;
% fname = "D:\2019.07.10\m190710.mat"; n=2;
% fname = "G:\2019.12.18\m191218.mat"; n=1;
% fname = "G:\2019.12.11\m191211.mat"; n=3;
fname = "H:\2020.11.30\m201130.mat"; n=2; 
% cd('C:\Users\Ori\Desktop\Ori\2nd degree\mtdt');
addpath("C:\Users\orica\Dropbox\fcns_and_decript");
addpath('C:\Users\orica\Dropbox\master degree\codes');
addpath("C:\Users\Ori\Dropbox\fcns_and_decript")
addpath("C:\Users\Ori\Dropbox\master degree\codes")
global bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0 signal2_smooth basis params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result = struct('TSCA',struct,'Tmax',struct,'AOF',struct,'Corr',struct,'GLM',struct,'Nadav',struct);
[result.TSCA.maps,result.Tmax.maps,result.AOF.maps,result.Corr.maps,result.GLM.maps,result.Nadav.maps,tscaBoth] = GenerateMaps(fname,n,9); % 3rd parameter: what, 8=loc8,9=loc9,92=mvngbars2hz
% <----- generate retinotopic maps from the individual maps
fn = fieldnames(result);
for i=1:length(fn) % iterate the methods 
    [~,result.(fn{i}).retinotopicMap] = retinotopicMapFromIndividualMaps(result.(fn{i}).maps,1,fn{i},80);
end
[~,resultBoth] = retinotopicMapFromIndividualMaps(tscaBoth.withGLM,1,'GLM+TSCA');
% ---- RUN UNTIL HERE -----
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
fname = "G:\181218\m181218.mat"; n=2;
% fname = "E:\200121\m200121.mat"; n=5;
[cf1 cfn trsh0]=strt_up(fname, n);  


result = Summary2.result;
params = Summary2.params;
% [cf1 cfn trsh0]=strt_up(fname, n);  
fn = fieldnames(result);
retmapsAll = figure("name","retmapsAll");
% indmapsAll = figure("name","indmapsAll");
Titles = {'(1)','(2)','(3)','(4)','(5)','(6)','(7)','(8)'};
for i=1:length(fn)-1 % iterate the methods    
    figure(retmapsAll); subplot(2,4,i)
    [~,r] = retinotopicMapFromIndividualMaps(result.(fn{i}).maps,0,fn{i},95);
    imf2(r);
%     imagesc(rshp(result.(fn{i}).retinotopicMap));
    title(Titles{i});
%     if contains(fn{i},'nadav','ignorecase',true)
%         title('M.P.T');
%     elseif contains(fn{i},'tscano','ignorecase',true)
%         title('TSCA');
%     elseif contains(fn{i},'tscaw','ignorecase',true)
%         title('TSCA & GLM');
%     else
%         title(fn{i});
%     end
%     r = plotMaps(result.(fn{i}).maps,fn{i},1);
end
subplot(2,4,8)
title(Titles{end});
% with and without GLM
% figure;
% subplot 121
% [~,r] = retinotopicMapFromIndividualMaps(result.TSCAnoGLM.maps,0,'TSCA without GLM',92);
% imf2(r); title('TSCA without GLM');
% subplot 122
% [~,r] = retinotopicMapFromIndividualMaps(result.TSCAwGLM.maps,0,'TSCA with GLM',92);
% imf2(r);title('TSCA after GLM denoising');
%% plot average response of 1 condition
t = linspace(0,1,100);
figure;
plot(t,params.experiment.responseSig)
ylabel('amp [au]');xlabel('time [sec]');
avgRes = zeros(1,params.experiment.T1);
zz = params.experiment.ZZ;
zz3d = rshp(zz);
avgAll = MinMaxNorm(mean(zz3d,3));
brn3d = reshape(repmat(MinMaxNorm(brn),1,1,3),[],3);
% avgAll(avgAll(:,1)<prctile(avgAll(:,1),95),:) = brn3d(avgAll(:,1)<prctile(avgAll(:,1),95),:);
figure;imf2(brn3d);hold on;imf2(rshp(avgAll),prctile(avgAll(:),93));
zzROI = [];
for i =1:params.experiment.N
    thisMap = result.TSCAwGLM.maps(:,:,i);
    [row,col] = find(thisMap>0.8);
    thisZZ = zz(:,(i-1)*params.experiment.T1+1:i*params.experiment.T1);
    zzROI = [zzROI,mean(thisZZ)];
    avgRes = avgRes + mean(thisZZ(sub2ind(size(brn),row,col),:));
end
figure;
% yyaxis left
plot(t,avgRes/params.experiment.N);ylabel('Z score');
yyaxis right
plot(t,params.experiment.responseSig)
ylabel('amp [au]');xlabel('time [sec]');
legend('Experimental Average Response','Theoretical Response');
figure;plot(0:0.01:8-0.01,zzROI);ylabel('Z score');xlabel('time [sec]');