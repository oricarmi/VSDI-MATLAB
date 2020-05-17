clear all; close all; clc;
% cd('C:\Users\Ori\Desktop\Ori\2nd degree\mtdt');
addpath("C:\Users\orica\Dropbox\fcns_and_decript");
addpath('C:\Users\orica\Dropbox\master degree\codes');
% fname='C:\Users\Ori\Desktop\Ori\2nd degree\mtdt\180905\m180905.mat';
% fname='C:\Users\Ori\Desktop\Ori\2nd degree\mtdt\191119\m191119.mat';
% fname='C:\Users\Ori\Desktop\Ori\2nd degree\mtdt\191212\m191212.mat';
% fname = "C:\Users\Ori\Desktop\Ori\2nd degree\mtdt\191119\m191119.mat"
fname = "G:\2020.01.21\m200121.mat";
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
n_strt=3; px=3; lgn00=lgn(n_strt:end,:);

% clz;
[rf cf]=avg_mx(cfn(n_strt:end), px); %find maximum of average of conditions, 2nd input is spatial filter width


%% pick area to work on from all conditions
% clz; 
cfn0=cfn; sz0=sz; brn0=brn;
cfn=pk_zn2(cfn,3);
%% average signal of all conditions
signal = cell2mat([cellfun(@mean,cfn(2:end),'UniformOutput',false)]');
signal2 = mean(signal);
signal2_smooth = filtfilt(ones(1,7),7,signal2); figure;plot(signal2_smooth);
% cfn=cfn0; sz=sz0; brn=brn0; %restore to original condition
%% tsca 
noises(1) = struct('f0',3.3,'bw',0.1); % HR?
noises(2) = struct('f0',22.85,'bw',0.1); % HR?
noises(3) = struct('f0',12.79,'bw',0.1); % Breathing? 
% noises(4) = struct('f0',0.68,'bw',0.1); % Breathing? 
% [status,img,fitted2DG,fitResult] = tsca_loc8(noises,[0.1 0.3],[2 10],[1 -1.*ones(1,length(noises)+1)],3,1,3,1,0);
[status,img,fitted2DG,fitResult] = tsca_loc8_2(noises,[],[2 10],[1 -1.*ones(1,length(noises)+1)],20,1,2.5,1,0);


%% calculation by percentile/threshold                           range&average
vc(2)=0; xs=[0.3 1]; %% remove limits and title from plots if vc(2)=1 and show scale xs=[x y]
vc(3)=1; %show reject/anti wave if vc(3)=1;
clz; plt_on=0; cnd=[]; p=1; x=1; t_lmts=[0.1 0.3];  stl=[3 50]; t_flt=10; s_flt=[10 2];  fgn=100; scl=[]; cmap=colormap(jet);
[mxc rxc]=xp_mp(cfn, cnd, p, x, t_lmts, stl,  t_flt, s_flt);
%figure(99) will give you the color index/scale

%% plot 8 or 4 location figures in right order/position
loc84pt(500, cfn, 0.1);

%%%red%%%%green%%blue%%magenta%cyan%%yellow%white%black

%%% pc=0 center by maximum amplitude, otherwise center by centroid  
clz; pc=0; [ds rf cf L]=cmbn84(mxc(:,2:end-1), rxc(:,2:end-1), 2); %% combine 8 or 4 conditions
