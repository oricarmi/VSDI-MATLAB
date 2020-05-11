%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERAL %%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc;
% fname = "G:\2020.01.21\m200121.mat"; n=6;
% fname = "G:\191119\m191119.mat"; n=2;
fname = "G:\180904\m180904.mat"; n=4;
% fname = "D:\2019.07.10\m190710.mat"; n=3;
% cd('C:\Users\Ori\Desktop\Ori\2nd degree\mtdt');
addpath("C:\Users\orica\Dropbox\fcns_and_decript");
addpath('C:\Users\orica\Dropbox\master degree\codes');
global bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0 signal2_smooth basis params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result = struct('TSCA',struct,'Tmax',struct,'AOF',struct,'Corr',struct,'GLM',struct,'Nadav',struct);
[result.TSCA.maps,result.Tmax.maps,result.AOF.maps,result.Corr.maps,result.GLM.maps,result.Nadav.maps] = GenerateMaps(fname,n,8); % 3rd parameter: what, 8=loc8,9=loc9,92=mvngbars2hz
% <----- generate retinotopic maps from the individual maps
fn = fieldnames(result);
for i=1:length(fn) % iterate the methods 
    [~,result.(fn{i}).retinotopicMap] = retinotopicMapFromIndividualMaps(result.(fn{i}).maps,8,fn(i));
end
% ---->
[result.TSCA.performance,result.Tmax.performance,result.AOF.performance,result.Corr.performance,result.GLM.performance,result.Nadav.performance] = performanceRealData(result);
Summary = struct('params',params,'result',result','description','loc 8 vis,191119 n=4'); 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% LM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = 0:1/fs:(1-1/fs);
alphas = [50 200; 10 100; 10 50;50 150]./1000;
randomAlphas = [...
    alphas(1,1)+rand(15000,1)*(alphas(1,2)-alphas(1,1))...
    alphas(2,1)+rand(15000,1)*(alphas(2,2)-alphas(2,1))...
    alphas(3,1)+rand(15000,1)*(alphas(3,2)-alphas(3,1))...
    alphas(4,1)+rand(15000,1)*(alphas(4,2)-alphas(4,1))];
r = zeros(size(randomAlphas,1),length(t));
for i=1:size(r,1) % iterate random signals
    for j=1:size(r,2) % iterate time
        if t(j)<=randomAlphas(i,1) || t(j)>=(randomAlphas(i,1)+randomAlphas(i,2)+randomAlphas(i,3)+randomAlphas(i,4)) % if before latency or after decay, it is zero
            continue
        elseif randomAlphas(i,1)<=t(j) && t(j)<=(randomAlphas(i,1)+randomAlphas(i,2)) % if in rise time
            r(i,j) = 0.5*(1-cos(pi*(t(j)-randomAlphas(i,1))/randomAlphas(i,2)));
        elseif (randomAlphas(i,1)+randomAlphas(i,2))<=t(j) && t(j)<=(randomAlphas(i,1)+randomAlphas(i,2)+randomAlphas(i,3)) % if in plateau time
            r(i,j) = 1;
        else %if (randomAlphas(i,1)+randomAlphas(i,2)+randomAlphas(i,3))<=t(j) && t(j)<=(randomAlphas(i,1)+randomAlphas(i,2)+randomAlphas(i,3)+randomAlphas(i,4)) % if in decaying time
            r(i,j) = 0.5*(1+cos(pi*(t(j)-randomAlphas(i,1)-randomAlphas(i,2)-randomAlphas(i,3))/randomAlphas(i,4)));
        end
    end
end
% figure;plot(r(randomAlphas(:,3)<0.03,:)');
[U,S,V] = svd(r,'econ');
basis = V(:,1:20);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% TSCA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [mapTSCA,fitted2DG,correlation] = VSDI_Analysis(fname,n,'MB2HZ',[],1); 
[mapTSCA,fitted2DG,correlation] = VSDI_Analysis(fname,n,8,[],1); 
retinotopicMapTSCA = retinotopicMapFromIndividualMaps(cat(3,mapTSCA{:}),[],2);
% for i =1:8
%     mapTSCA2{i} = MinMaxNorm(mapTSCA{i});
% end
% retinotopicMapTSCA = retinotopicMapFromIndividualMaps(cat(3,mapTSCA2{:}),2);
% retinotopicMapFitted = retinotopicMapFromIndividualMaps(cat(3,fitted2DG{:}),2);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%% Nadav analysis %%%%%%%%%%%%%%%%%%%%%%%%%%
%% upload data
clz;  [cf1 cfn trsh0]=strt_up(fname, n);    
Z = []; % concatenate
for i=4:length(cfn)
    Z = [Z cfn{i}];
end

%% pick area to work on from all conditions
%clz; 
cfn0=cfn; sz0=sz; brn0=brn; % run only on first time
cfn=pk_zn2(cfn,1);
%% params
% cfn=cfn0; sz=sz0; brn=brn0; %restore to original condition
% calculation by percentile/threshold                           range&average
vc(2)=0; xs=[0.3 1]; %% remove limits and title from plots if vc(2)=1 and show scale xs=[x y]
vc(3)=1; %show reject/anti wave if vc(3)=1;
%settle within stl(1) std at stl(2) pixels from the end
% s_flt(1)= gaussian sigman and s_flt(2) percent cutoff
% clz; plt_on=1; cnd=[]; p=[]; x=0; t_lmts=[0 1];  stl=[3 10]; t_flt=0; s_flt=[5 10];  fgn=100; scl=[]; cmap=colormap(jet);
clz; plt_on=1; cnd=4:12; p=[]; x=0; t_lmts=[0 1];  stl=[1 10]; t_flt=0.5; s_flt=[5 10];  fgn=100; scl=[]; cmap=colormap(jet);
[mxc rxc]=xp_mp(cfn, cnd, p, x, t_lmts, stl,  t_flt, s_flt); %figure(99) will give you the color index/scale
loc9p(500, cfn(3:end), 0.1); %% loc9 plot
clz; pc=0; [ry cx mx]=cmbn9(mxc(:,2:end)); %% combine 9 conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%% Average of frames %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = []; % concatenate
for i=1:length(cfn)
    if i==2 
        continue
    end
    Z = [Z cfn{i}];
end
ZZ = preProcess(reshape(Z,[],800),800);
ZZ = rshp(Z);
for i=1:8
%     mapAOF{i} = imgaussfilt(mean(ZZ(:,:,(i-1)*50+6:(i-1)*50+20),3),3);
    mapAOF{i} = imgaussfilt(mean(ZZ(:,:,(i-1)*100+6:(i-1)*100+20),3),0.8);
end
plotMaps(cat(3,mapAOF{:}),8,0);
retinotopicMapAOF = retinotopicMapFromIndividualMaps(cat(3,mapAOF{:}),2);
%% %%%%%%%%%%%%%%%%%%%%%%%%%% T max %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Z = []; % concatenate
for i=3:length(cfn)
    Z = [Z cfn{i}];
end
ZZ = preProcess(Z,size(Z,2));
mapT = Tmax(reshape(ZZ,[],800),8,0,0.5); % signal, what (8,9,anything) , std for r, std for imgaussfilt
retinotopicMapTmax = retinotopicMapFromIndividualMaps(cat(3,mapT),2);
% figure;imagesc(ump.*[0:size(brn,2)-1]./1000, ump.*[0:size(brn,1)-1]./1000,reshape(tmax2,270,[],3));colormap(jet); 
% figure;imagesc(reshape(colors,8,1,3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%% Corr %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mapCorr = Correlation2All(reshape(ZZ,[],800),8,[]);
retinotopicMapTmax = retinotopicMapFromIndividualMaps(mapCorr,2);

