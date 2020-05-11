clear all; close all; clc;
cd('C:\Users\Ori\Desktop\Ori\2nd degree\mtdt');
addpath('C:\Users\Ori\Dropbox\fcns_and_decript');
addpath('C:\Users\Ori\Dropbox\master degree\codes');
% fname='C:\Users\Ori\Desktop\Ori\2nd degree\mtdt\190221\m190821.mat';
fname = "E:\2020.01.22\m200122.mat";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% copy&paste below%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0 signal2_smooth

% vc=vector of conditions:
%1) loc428/strt_up/loc84pt control 8loc vs 42loc
%2) remove title/limits from xp_mp
%3) plot antiwave in xp_mp

% cfn=cfn0; sz=sz0; brn=brn0; restore to original condition

%% upload data
clz; n=4; [cf1 cfn trsh0]=strt_up(fname, n);
%% run till here

%% initial nlz
% last input: df/f=0, std=1
[dt_dfof dt_nrm]=int_nlz(cf1, cfn, n, 0.1, 1);

%% show avg
n_strt=2; px=3; lgn00=lgn(n_strt:end,:);

clz; [rf cf]=avg_mx(cfn(n_strt:end), px); %find maximum of average of conditions, 2nd input is spatial filter width


%% pick area in figure to show waveforms:

fgr_n=1 %% fgr_n=figure number

i_sub=pk_r(cfn, cf1, fgr_n);


%% pick area to work on from all conditions
%clz; 
cfn0=cfn; sz0=sz; brn0=brn; % run only on first time
cfn=pk_zn2(cfn,1);
% cfn=cfn0; sz=sz0; brn=brn0; %restore to original condition

%% LOC9 NIR - No parameters
% calculation by percentile/threshold                           range&average
vc(2)=0; xs=[0.3 1]; %% remove limits and title from plots if vc(2)=1 and show scale xs=[x y]
vc(3)=1; %show reject/anti wave if vc(3)=1;
%settle within stl(1) std at stl(2) pixels from the end
% s_flt(1)= gaussian sigman and s_flt(2) percent cutoff
clz; plt_on=1; cnd=[]; p=[]; x=0; t_lmts=[0 1];  stl=[3 10]; t_flt=0; s_flt=[5 10];  fgn=100; scl=[]; cmap=colormap(jet);
[mxc rxc]=xp_mp(cfn, cnd, p, x, t_lmts, stl,  t_flt, s_flt); %figure(99) will give you the color index/scale
loc9p(500, cfn, 0.1); %% loc9 plot

clz; pc=0; [ry cx mx]=cmbn9(mxc(:,2:end)); %% combine 9 conditions

%% LOC9 NIR
% calculation by percentile/threshold                           range&average
vc(2)=0; xs=[0.3 1]; %% remove limits and title from plots if vc(2)=1 and show scale xs=[x y]
vc(3)=1; %show reject/anti wave if vc(3)=1;
%settle within stl(1) std at stl(2) pixels from the end
% s_flt(1)= gaussian sigman and s_flt(2) percent cutoff
% clz; plt_on=1; cnd=[]; p=[]; x=0; t_lmts=[0 1];  stl=[3 10]; t_flt=0; s_flt=[5 10];  fgn=100; scl=[]; cmap=colormap(jet);
clz; plt_on=1; cnd=[]; p=[]; x=0; t_lmts=[0 1];  stl=[1 10]; t_flt=0.5; s_flt=[5 10];  fgn=100; scl=[]; cmap=colormap(jet);
[mxc rxc]=xp_mp(cfn, cnd, p, x, t_lmts, stl,  t_flt, s_flt); %figure(99) will give you the color index/scale
loc9p(500, cfn(3:end), 0.1); %% loc9 plot

clz; pc=0; [ry cx mx]=cmbn9(mxc(:,2:end)); %% combine 9 conditions


c=cnv2cll(mxc);
%B C DL UR D U DR UL R L
%1 2 3  4  5 6 7  8  9 10
u=c{8}+c{6}+c{4};
d=c{3}+c{5}+c{7};
r=c{4}+c{7}+c{9};
l=c{3}+c{8}+c{10};
ul=c{6}+c{8}+c{10};
ur=c{6}+c{8}+c{10};
dl=c{10}+c{3}+c{5};
dr=c{9}+c{5}+c{7};
s=c{3}+c{4}+c{5}+c{6}+c{7}+c{8}+c{9}+c{10};
cfnl=l-r; iout=imf(cfnl); colormap gray; title('L-R') %figure; hist(iout(:));
cfnl=u-d; iout=imf(cfnl); colormap gray; title('U-D')%figure; hist(iout(:));
cfnl=ul-dr; iout=imf(cfnl); colormap gray; title('UL-DR') %figure; hist(iout(:));
cfnl=ur-dl; iout=imf(cfnl); colormap gray; title('UR-DL')%figure; hist(iout(:));
cfnl=c{2}-s; iout=imf(cfnl); colormap gray; title('C-S')%figure; hist(iout(:));
%% weighted average
cfnl=(l-r)./(l+r); iout=imf(cfnl); colormap gray; title('L-R/L+R') %figure; hist(iout(:));
cfnl=(u-d)./(u+d); iout=imf(cfnl); colormap gray; title('U-D/U+D')%figure; hist(iout(:));
cfnl=(ul-dr)./(ul+dr); iout=imf(cfnl); colormap gray; title('UL-DR/UL+DR') %figure; hist(iout(:));
cfnl=(ur-dl)./(ur+dl); iout=imf(cfnl); colormap gray; title('UR-DL/UR+DL')%figure; hist(iout(:));
cfnl=(c{2}-s)./(c{2}+s); iout=imf(cfnl); colormap gray; title('C-S/C+S')%figure; hist(iout(:));
%% 
%% average signal of all conditions
signal = cell2mat([cellfun(@mean,cfn(2:end),'UniformOutput',false)]');
signal2 = mean(signal);
signal2_smooth = filtfilt(ones(1,7),7,signal2); figure;plot(signal2_smooth);
% cfn=cfn0; sz=sz0; brn=brn0; %restore to original condition
%% tsca 
noises(1) = struct('f0',3.125,'bw',0.1); % HR?
% noises(2) = struct('f0',3.125,'bw',0.2); % HR?
% noises(2) = struct('f0',0.68,'bw',1); % Breathing? 
[status,img,fitted2DG,fitResult] = tsca_loc9(noises,[0.1],[2 10],[1 -1.*ones(1,length(noises)+1)],3,1,3,1,0);
