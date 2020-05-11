function [img,fitted2DG,correlatoin,fitResult,offsets] = VSDI_Analysis(fname,n,what,tlim,toFit,all)
% cd('C:\Users\orica\Desktop\Ori\2nd degree\mtdt');
addpath('C:\Users\orica\Dropbox\fcns_and_decript');
addpath('C:\Users\orica\Dropbox\master degree\codes');
% fname='C:\Users\Ori\Desktop\Ori\2nd degree\mtdt\180905\m180905.mat';
% fname='C:\Users\Ori\Desktop\Ori\2nd degree\mtdt\191119\m191119.mat';
% fname='C:\Users\Ori\Desktop\Ori\2nd degree\mtdt\191212\m191212.mat';
% fname = 'F:\2019.12.18\m191218';
% fname = "F:\2019.12.18\MAT\E0\NI232_E0B0.BLK.mat";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% copy&paste below%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<3
    disp('must provide fname and n and what');
    return
elseif nargin<4
    tlim = [];
    all = 0;
    toFit = 0;
elseif nargin<5
    all = 0;
    toFit = 0;
elseif nargin<6
    all = 0;
end
global bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0 signal2_smooth basis
clz; [cf1 cfn trsh0]=strt_up(fname, n);   

% <---- select area of response to build response signal
n_strt=2; px=3; lgn00=lgn(n_strt:end,:);
[rf cf]=avg_mx(cfn(n_strt:end), px); %find maximum of average of conditions, 2nd input is spatial filter width
cfn0=cfn; sz0=sz; brn0=brn;
% cfn=pk_zn2(cfn,3);
% signal = cell2mat([cellfun(@mean,cfn(2:end),'UniformOutput',false)]');
% signal2 = mean(signal); x = [1:length(signal2)];
% signal2 = filtfilt(ones(1,7),7,signal2); %signal2(40:end) = 0;
% figure;hold on;plot(x,signal2);
% ver2use = input('which version to use? [1,2]');
% if ver2use==1
%     gaussEqn = 'a1*exp(-((x-b1)/c1)^2)+a2*exp(-((x-b2)/c2)^2)';
%     fittedExp = fit(x',signal2',gaussEqn,'start',[max(signal2)/2 round(length(signal2)/3) 5 max(signal2) round(length(signal2)/2) 14]);
%     signal2_smooth = fittedExp.a1*exp(-((x-fittedExp.b1)/fittedExp.c1).^2)+fittedExp.a2*exp(-((x-fittedExp.b2)/fittedExp.c2).^2);
% %     signal2_smooth = polyfit(x',signal2',6);% fit(x',signal2',gaussEqn,'exclude',x>30);
% %     signal2_smooth = polyval(signal2_smooth,x);
% else
%     signal2_smooth = zeros(1,length(x));
%     [~,ind] = max(diff(signal2));
%     signal2_smooth(ind:ind+5) = 1;
% end
% % signal2_smooth = signal2;
% plot(x,signal2_smooth);
% % ---->
% < ---- pick area of blood vessels to determine noises
% cfn=cfn0; sz=sz0; brn=brn0; %restore to original condition
cfn=pk_zn2(cfn,2);
Z = []; % concatenate
for i=1:length(cfn)
    if i==2 && any(what==8)
        continue;
    end
    Z = [Z cfn{i}];
end
Z = Z - mean(Z,2);

% [ppx,f] = periodogram(mean(Z),[],[],fs);
[ppx,f] = pmtm(mean(Z),2,pow2(nextpow2(length(mean(Z)))),fs); % calculate spectrum
figure; % show time course and spectrum
subplot(1,2,1);plot(mean(Z)); xlabel('frame'); ylabel('amp'); title('time course of average signal');
subplot(1,2,2); plot(f,ppx); xlabel('frequency [Hz]');ylabel('amp');title('spectrum of average signal');
prompt = 'What frequencies are noise?[frq1,frq2,...]';
f0 = input(prompt);
if length(f0)>0
    for i=1:length(f0)
        noises(i) = struct('f0',f0(i),'bw',0.1); 
    end
else
    noises = [];
end
% ---->
cfn=cfn0; sz=sz0; brn=brn0; %restore to original condition
prompt = 'Whole Brain? [0,1]';
allBrn = input(prompt);
if allBrn == 0
   cfn=pk_zn2(cfn,2);
end
% < --- concatenate all 10 seconds to run in video
Z = [];
for i=1:length(cfn)
    if i==2 && any(what==8)
        continue;
    end
    Z = [Z cfn{i}];
end
ZZZ = preProcess(Z,size(Z,2));
ZZZ = reshape(ZZZ,size(Z,1),[]);
% ZZZ = Z - min(Z(:));
implay(rshp(Z),20); % play video to check results of TSCA 
% ---->
% < --- run TSCA
switch what
    case 8
        [~,img,fitted2DG,correlatoin,fitResult,offsets] = tsca_loc8_2(ZZZ(:,101:end),noises,tlim,[2 10],[1 -0.05.*ones(1,length(noises)+1)],49,1,0.8,toFit,0);
    case 9
        [~,img,fitted2DG,correlatoin,fitResult,offsets] = tsca_loc9_2(ZZZ(:,101:end),noises,tlim,[2 10],[1 -0.05*ones(1,length(noises)+1)],49,1,1.5,toFit,0);
    case 'MB2HZ'
        [~,img,fitted2DG,correlatoin,fitResult,offsets] = tsca_movingBars_2hz(ZZZ(:,151:end),noises,tlim,[2 10],[1 -0.01*ones(1,length(noises)+1)],49,1,1.5,toFit,0);
    otherwise 
        [~,img,fitted2DG,correlatoin,fitResult] = tsca_ff(noises,tlim,[2 10],[1 -0.01.*ones(1,length(noises)+1)],20,1,0,0);
end
% ---->
% <---- NADAV analysis
if all
    [dt_dfof dt_nrm]=int_nlz(cf1, cfn, n, 0.1, 1);
    cfn=pk_zn2(cfn,3);
    vc(2)=0; xs=[0.3 1]; %% remove limits and title from plots if vc(2)=1 and show scale xs=[x y]
    vc(3)=1; %show reject/anti wave if vc(3)=1;
    plt_on=1; cnd=[]; p=1; x=1; t_lmts=[0.1 0.3];  stl=[3 50]; t_flt=10; s_flt=[10 2];  fgn=100; scl=[]; cmap=colormap(jet);
    [mxc rxc]=xp_mp(cfn, cnd, p, x, t_lmts, stl,  t_flt, s_flt);
    loc84pt(500, cfn, 0.1);
    pc=0; [ds rf cf L]=cmbn84(mxc(:,2:end-1), rxc(:,2:end-1), 2); %% combine 8 or 4 conditions
end
% ---->
end

