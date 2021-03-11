clear all; close all; clc;
% fname = "C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\2021.01.18\m210118.mat"; % path to .mat file (mYYMMDD)
fname="D:\dataForComparison\181218\m181218.mat"
% fname ="H:\2021.02.09\m210209.mat";

% fname ="C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\2021.01.18 - target flankers\m210118.mat";
n=2; % experiment number 
what = 8; % type of experiment. 8 - loc 8. 9 - loc 9 grid. 92 - 9 moving bars/diagonal dots, 2Hz
% cd('C:\Users\Ori\Desktop\Ori\2nd degree\mtdt');
addpath("C:\Users\orica\Dropbox\fcns_and_decript");
addpath('C:\Users\orica\Dropbox\master degree\codes');
addpath("C:\Users\Ori\Dropbox\fcns_and_decript")
addpath("C:\Users\Ori\Dropbox\master degree\codes")
global bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0 signal2_smooth basis params
[cf1 cfn trsh0]=strt_up(fname, n); % parse .mat file (mYYMMDD) (legacy code)
load('C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\comparision results 2\loc 8,181218 n=2.mat')
%%
t = linspace(0,size(Z,2)/100,size(Z,2));
mapTSCA = Summary2.result.TSCAwGLM.maps;
ZZ = Summary2.params.experiment.ZZ;
avgTargetSpatial = mean(mapTSCA,3);
figure;imf2(rshp(tempbrn)); hold on; imf2(avgTargetSpatial,0.4);
[row,col] = find( avgTargetSpatial>=0.4);
targetROI = zeros(sz); 
for i = 1:length(row)
    targetROI(row(i),col(i)) = 1; 
end
figure;
subplot 121
imf2(targetROI);title('target ROI');
indTargetROI = sub2ind(sz,row,col);
subplot 122
plot(t,mean(ZZ(indTargetROI,:));