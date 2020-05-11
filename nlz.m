clear all; close all; clc;
cd('C:\Users\Ori\Desktop\Ori\2nd degree\mtdt');
addpath('C:\Users\Ori\Dropbox\fcns_and_decript');
% fname='C:\Users\Ori\Desktop\Ori\2nd degree\mtdt\180905\m180905.mat';
addpath('C:\Users\Ori\Dropbox\master degree\codes');
fname='C:\Users\Ori\Desktop\Ori\2nd degree\mtdt\191212\m191212.mat';
% fname = 'F:\2019.12.12\m191212.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% copy&paste below%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global bsln fs sz ump rot fgn brn frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms

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

%% combine experiment 2 loc4 to 1 loc8 experiment
[cfn n]=cmbn428(c_f, 2 , 3);

%% localized waveforms    (several maps with local response - in a matrix)
clz; wvs=wv_loc(cfn, 6);

%% raw maps 36 is number of frames

%scl=0 local scaling


%scl=1 global (matrix) scaling


%scl=2 global (cell) scaling


%scl=[mn mx] choose scaling

%%%%%%%%%%%%%%%%%%cplt(m0, v, scl, rc)
n_frame=9+[1:9]

clz;  cplt(cfn, n_frame, 0); % this enable to make maps with various scaling - IMPORTANT!!

%%combine raw maps (color)

%% red=1    green=2    blue=3   magenta=4   cyan=5 yellow=6    white=7   black =8
% n+1 color is the combination of multiple conditions
x=0.1; %% x= threshold

n_frame=7+[1:9]; %% frames to show

clz; im_out=ccmb(cfn(3:end), n_frame, x); %%% n_frame number of frames


%% pick area to work on from all conditions
clz; cfn0=cfn; sz0=sz; brn0=brn;
cfn=pk_zn(cfn);

% cfn=cfn0; sz=sz0; brn=brn0; %restore to original condition


%% calculation by percentile/threshold                           range&average
vc(2)=0; xs=[0.3 1]; %% remove limits and title from plots if vc(2)=1 and show scale xs=[x y]
vc(3)=1; %show reject/anti wave if vc(3)=1;
clz; plt_on=1; cnd=[]; p=1; x=0; t_lmts=[0 1];  stl=[3 25]; t_flt=10; s_flt=[10 5];  fgn=100; scl=[]; cmap=colormap(jet);
[mxc rxc]=xp_mp(cfn, cnd, p, x, t_lmts, stl,  t_flt, s_flt);
%figure(99) will give you the color index/scale

%% open/no prms
% calculation by percentile/threshold                           range&average
vc(2)=0; xs=[0.3 1]; %% remove limits and title from plots if vc(2)=1 and show scale xs=[x y]
vc(3)=1; %show reject/anti wave if vc(3)=1;
%settle within stl(1) std at stl(2) pixels from the end
% s_flt(1)= gaussian sigman and s_flt(2) percent cutoff
clz; plt_on=1; cnd=[]; p=[]; x=0; t_lmts=[0 1];  stl=[100 10]; t_flt=0; s_flt=[0 3];  fgn=100; scl=[]; cmap=colormap(jet);
[mxc rxc]=xp_mp(cfn, cnd, p, x, t_lmts, stl,  t_flt, s_flt);
%figure(99) will give you the color index/scale

%calculation by percentile/threshold                           range&average
vc(2)=1; xs=[0.3 1]; %% remove limits and title from plots if vc(2)=1 and show scale xs=[x y]
vc(3)=1; %show reject/anti wave if vc(3)=1;
clz; plt_on=1; cnd=[]; p=0.1; x=4; t_lmts=[0.15 0.35];  stl=[1 10]; t_flt=10; s_flt=[0 5];  fgn=100; scl=[]; cmap=colormap(jet);
[mxc rxc]=xp_mp(cfn, cnd, p, x, t_lmts, stl,  t_flt, s_flt);
%figure(99) will give you the color index/scale

%% NIR prms
%calculation by percentile/threshold                           range&average
%settle within stl(1) std at stl(2) pixels from the end
% s_flt(1)= gaussian sigman and s_flt(2) percent cutoff
vc(2)=0; xs=[0.3 1]; %% remove limits and title from plots if vc(2)=1 and show scale xs=[x y]
vc(3)=1; %show reject/anti wave if vc(3)=1;
clz; plt_on=1; cnd=[]; p=[]; x=0; t_lmts=[0.1 0.2];  stl=[1 10]; t_flt=10; s_flt=[5 10];  fgn=100; scl=[]; cmap=colormap(jet);
[mxc rxc]=xp_mp(cfn, cnd, p, x, t_lmts, stl,  t_flt, s_flt);
%figure(99) will give you the color index/scale

%% plot 8 or 4 location figures in right order/position
loc84pt(500, cfn, 0.1);

%% VIR prms
%calculation by percentile/threshold                           range&average
%settle within stl(1) std at stl(2) pixels from the end
% s_flt(1)= gaussian sigman and s_flt(2) percent cutoff
vc(2)=0; xs=[0.3 1]; %% remove limits and title from plots if vc(2)=1 and show scale xs=[x y]
vc(3)=0; %show reject/anti wave if vc(3)=1;
clz; plt_on=1; cnd=[]; p=1; x=0; t_lmts=[0.01 0.5];  stl=[3 10]; t_flt=10; s_flt=[3 10];  fgn=100; scl=[]; cmap=colormap(jet);
[mxc rxc]=xp_mp(cfn, cnd, p, x, t_lmts, stl,  t_flt, s_flt);
%figure(99) will give you the color index/scale

%% LOC9 NIR
% calculation by percentile/threshold                           range&average
vc(2)=0; xs=[0.3 1]; %% remove limits and title from plots if vc(2)=1 and show scale xs=[x y]
vc(3)=1; %show reject/anti wave if vc(3)=1;
%settle within stl(1) std at stl(2) pixels from the end
% s_flt(1)= gaussian sigman and s_flt(2) percent cutoff
clz; plt_on=1; cnd=[]; p=[]; x=0.5; t_lmts=[0.1 0.5];  stl=[1 10]; t_flt=0; s_flt=[5 10];  fgn=100; scl=[]; cmap=colormap(jet);
[mxc rxc]=xp_mp(cfn, cnd, p, x, t_lmts, stl,  t_flt, s_flt); %figure(99) will give you the color index/scale
loc9p(500, cfn, 0.1); %% loc9 plot

clz; pc=0; [ry cx mx]=cmbn9(mxc(:,2:end)); %% combine 9 conditions





%% plot 8 or 4 location figures in right order/position
loc84pt(500, cfn, 0.1);

%%%red%%%%green%%blue%%magenta%cyan%%yellow%white%black

%%% pc=0 center by maximum amplitude, otherwise center by centroid  
clz; pc=0; [ds rf cf L]=cmbn84(mxc(:,2:end-1), rxc(:,2:end-1), 2); %% combine 8 or 4 conditions

%% combine reference images
 fg1='D:\2018.08.28\180828\E3pre_1.bmp';
 fg2='D:\2018.08.28\180828\E3pre_2.bmp';
 
 fg1='D:\2018.09.05\MAT\E3post_FF.bmp';
 fg2='D:\2018.09.05\MAT\E3post_loc8.bmp';
 
 fg1='D:\2018.08.22\MAT\E0pre.bmp';
 fg2='D:\2018.08.22\MAT\E0post.bmp';
 
 fg1='D:\2018.09.05\MAT\E2pre_LR.bmp';
 fg2='D:\2018.09.05\MAT\E2pre_RL.bmp';
 
  
 fg1="F:\2019.12.18\MAT\E1post_FF.BMP";
 fg2="F:\2019.12.18\MAT\E1pst_LOC8.BMP";


Cout=icmb(fg1, fg2);

%% stich experiments
fgnm=['fg180605_8.fig']; rcf=[rf{1} cf{1}]; smbl='3';
stch(rcf, smbl, fgnm);

%% degree/mm, retinal coordinates, cortcal COA
dg=15; crd=[0 0; 2 2; 0 4; 1 1; 2 4; 0 2; 2 0; 1 3]; in=1000.*ds{1}; %% ds=cortical map is in milimeters, in is in micrometers

dg=1; crd=[-2.349 1.83; -1.9 2.201; -2.391 2.538; -2.138 1.991; -1.928 2.566; -2.377 2.152; -1.893 1.823; -2.145 2.370];  in=1000.*ds{1}; %% ds=cortical map is in milimeters, in is in micrometers

dg=15; crd=[0 0; 0.5 0.5; 0.5 0; 0 0.5]; in=1000.*ds{1}; %% ds=cortical map is in milimeters, in is in micrometers

[cout rds]=cmf84(in,[rf{1} cf{1}], crd, dg); %% cmf in micrometer/degrees
%rds=[diagonal vertical horizontal] distance between centers


filename = '181114E0CMF.xlsx';
%A = [12.7 5.02 -98 63.9 0 -.2 56];
xlswrite(filename,cout)
%% Gaussian fit

mxn=0.1; %% mxn=lower threshold (any value below shows brain background)

clz; fgn=200; lgn00=lgn(3:end,:);

[df prm ds]=g2d(mxc(:,3:end), mxn);


%% pick area in figure to show waveforms:

fgr_n=1 %% fgr_n=figure number

i_sub=pk_r(cfn, cf1, fgr_n);

%% pick area around max point of picture to show waveforms

n_strt=2; px=3; lgn00=lgn(n_strt:end,:);

clz; [rf cf]=avg_mx(cfn(n_strt:end), px); %find maximum of average of conditions, 2nd input is spatial filter width

n_strt=3; lgn00=lgn(n_strt:end,:); px=3;
clz; [v_cntr vmx]=pk_mx([cf rf], cfn(n_strt:end), cf1(n_strt:end), px, 0); %% use avg_mx coordinates
Md=dmnx(v_cntr, [0 0.2]) %% calc mnmx difference
figure; [v_cntr vmx]=pk_mx([3.66 3.44], cfn(n_strt:end), cf1(n_strt:end), px, 0); %%pick spot by coordinates

nv=[1 3 5 7 9 10]; lgn00=lgn(nv,:); px=1;
figure; [v_cntr vmx]=pk_mx([2.44 3.36], cfn(nv), cf1(nv), px, 0); %%pick spot by coordinates
nv=[1 3 5 7 9 10]; lgn00=lgn(nv,:); px=5;
figure; [v_cntr vmx]=pk_mx([cf rf], cfn(nv), cf1(nv), px, 0); %%pick spot by coordinates

nv=[3 5 7 9]; lgn00=lgn(nv,:); px=1;
figure; [v_cntr vmx]=pk_mx([2.44 3.36], cfn(nv), cf1(nv), px, 10); %%pick spot by coordinates
figure; [v_cntr vmx]=pk_mx([cf rf], cfn(nv), cf1(nv), px, 0); %%pick spot by coordinates

nv=[1 3 5 7 9]; lgn00=lgn(nv,:); px=1;
figure; [v_cntr vmx bs]=pk_mx_bs([2.44 3.36], cfn(nv), cf1(nv), px, 0); %%pick spot by coordinates
figure; [v_cntr vmx bs]=pk_mx_bs([cf rf], cfn(nv), cf1(nv), px, 0); %%pick spot by coordinates




%% mathematical operations on data
%% choose mxc,  rxc, mlc or xtc
c=cnv2cll(mxc);
c34=c{3}-c{4};
c56=c{5}-c{6};
c36=c{3}-c{6};
figure; [ds rc]=rtmx([c34(:) c56(:) c36(:)]); %choose pixel with highest value for each condition

for k0=1:length(c); cn{k0}= c{k0}./mxa(c{k0}); end

cfnl=c{4}-c{3}; figure; iout=imbx(cfnl); colormap gray; figure; hist(iout(:));

cfnl=c{4}-c{3}; figure; imsp(cfnl,0.99); colormap gray

cfnl=cn{4}-(cn{3}+cn{5}+cn{6}); figure; imsp(cfnl,0.95); colormap gray

cfnl=cn{4}./cn{3}; figure; g=5; imsc(cfnl,g); colormap gray

cfnl=cn{4}./(cn{3}+cn{5}+cn{6}); figure; g=5; imsc(cfnl,g); colormap gray

cfnl=cn{3}-(cn{4}+cn{5}+cn{6}); figure; g=5; imsc(cfnl,g); colormap gray

cfnl=c{3}./(c{4}+c{5}+c{6}); figure; g=5; imsc(cfnl,g); colormap gray

cfnl=c{3}./(c{3}+c{4}+c{5}+c{6}); figure; g=5; imsc(cfnl,g); colormap gray

cfnl=c{2}./(c{2}+c{3}); figure; g=5; imsc(cfnl,g); colormap gray
cfnl=c{2}./mxa(c{2}+c{3}); figure; g=5; imsc(cfnl,g); colormap gray

cfnl=(c{3}+c{4}+c{5})./(c{6}+c{7}+c{8});

cfnl=c{3}-(c{4}+c{5}+c{6}+c{7}+c{8}+c{9}+c{10});

cfnl=c{4}-(c{3}+c{5}+c{6}+c{7}+c{8}+c{9}+c{10});

cfnl=c{5}-(c{3}+c{4}+c{6}+c{7}+c{8}+c{9}+c{10});

cfnl=1./((c{9}+c{4}+c{7})./(c{3}+c{8}+c{5}))

cfnl=(c{3}+c{9}+c{6})-(c{4}+c{8}+c{10}+c{7}+c{5})

cfnl=(c{10}+c{7}+c{5})-(c{3}+c{9}+c{6}+c{4}+c{8})
cfnl=(c{3}+c{5}+c{7}+c{9})-(4}+c{6}+c{8}+c{10})
cfnl=(c{3}-c{2})./(c{2}+c{3})
cfnl=c{3}./c{2};
cfnl=(c{3}./c{2})./(c{2}+c{3});
cfnl=c{5}./(c{3}+c{4}+c{5});

cfnl=c{7}-c{6};


%% scale figures aroung Gaussian mean with standard dev +-g
figure; g=1; imsc(cfnl,g); colormap gray
figure; img(cfnl); colormap gray

%% perfprm basic analysis on matrix of data
cc=cfn{2}-cfn{3};
clz; bsc_nlz(cc, 0.1);
%% raw maps 36 is number of frames

%scl=0 local scaling


%scl=1 global (matrix) scaling


%scl=2 global (cell) scaling


%scl=[mn mx] choose scaling

clz;  cplt({cc}, [8:16], 0);