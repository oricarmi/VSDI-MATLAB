%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath("C:\Users\orica\Dropbox\fcns_and_decript");
addpath('C:\Users\orica\Dropbox\master degree\codes');
global fs 
T = 1000; fs = 100;
P = 1600;
m = sqrt(P);
signals = struct;
noises = struct;
colors = [0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1];
stim_time = 10:20:90;
noiseSig = [0 0.1 0.3 0.5 1];
for kk = 1:length(noiseSig) % iterate different noise sig
runSummary = cell(4,20);
for k=1:20
%% construct signals
[I,J] = ndgrid(1:m,1:m); r = 4; 
% locs = [0.25 0.25; 0.25 0.75; 0.5 0.5; 0.75 0.25; 0.75 0.75]; 
locs = [0.4 0.4; 0.4 0.6; 0.5 0.5; 0.6 0.4; 0.6 0.6]; 
DC = 0.02; 
% figure;ind2plot = [1,5,9,13,17];
for i=1:5
    signals(i).time = 1*normpdf(0:0.1:(T-1)/10,stim_time(i),1);
    signals(i).space = imgaussfilt(double((I-m*locs(i,1)).^2+(J-m*locs(i,2)).^2<r^2),sqrt(r))+DC; % what is larger than r^2 is 1 (white), what is smaller is 0 (black)
%     subplot(3,6,ind2plot(i)) 
%     imagesc(signals(i).space); colormap(gray);
%     title(['sig. ' num2str(i) ' - spatial']);
%     subplot(3,6,ind2plot(i)+1)
%     plot(signals(i).time);
%     title(['sig. ' num2str(i) ' - temporal']);
%     xticks([0 100:200:T]);
end
%% construct noise
t = linspace(0,(T-1)./fs,T);
[I1,I2] = ndgrid([repmat(linspace(0,2*pi,m/2),1,2)]',[repmat(linspace(0,2*pi,m/2),1,2)]);
noises(1).time = normrnd(0,noiseSig(kk),T,1);
noises(1).space = cos(I1);
noises(2).time = noiseSig(kk)*cos(2*pi*1.*t'+rand(size(t')).*(2*pi));%normrnd(0,noiseSig,T,1);
noises(2).space = normrnd(0,1,m,m);
noises(3).time = noiseSig(kk)*cos(2*pi*3.*t'+rand(size(t')).*(2*pi));%normrnd(0,noiseSig,T,1);
noises(3).space = cos(I2);
% figure; ind2plot = [1,3,5];
% for i=1:length(noises)
%     subplot(3,2,ind2plot(i))
%     imagesc(noises(i).space); colormap(gray);
%     title(['noise ' num2str(i) ' - spatial']);
%     subplot(3,2,ind2plot(i)+1)
%     plot(noises(i).time);
%     title(['noise ' num2str(i) ' - temporal']);
% end
%% construct Z and show 9 frames
Z = zeros(m*m,T); % preallocate memory
for i = 1:T
    Z(:,i) = reshape( ...
        signals(1).time(i)*signals(1).space+signals(2).time(i)*signals(2).space+...
        signals(3).time(i)*signals(3).space+signals(4).time(i)*signals(4).space+...
        signals(5).time(i)*signals(5).space+...
        noises(1).time(i)*noises(1).space+...
        noises(2).time(i)*noises(2).space+...
        noises(3).time(i)*noises(3).space...
        ,[],1);
end
% figure; suptitle('the recorded signal (Z) at random frames');
% ind2show = sort([100;502;899;randi(1000,7,1)]);
% for i=1:length(ind2show)
%     subplot(2,5,i)
%     imagesc(reshape(Z(:,ind2show(i)),m,[])); colormap(gray);
%     title(['Frame #: ' num2str(ind2show(i))]);
% end
%% calc performance measures between theoretical signal and Z 
for i=1:length(signals)
    [~,I] = max(signals(i).time);
    mapORIG(:,:,i) = reshape(mean(Z(:,I-25:I+25),2),40,40);
    [origMSE(i),origPSNR(i),origCNR(i),origMSSIM(i),origCorr(i),origCP(i)] = getPerformance( mapORIG(:,:,i),signals(i).space,signals(i).space~=DC,signals(i).space==DC);
%     figure;imagesc(reshape(mean(Z(:,I-25:I+25),2),40,40));
end
retinotopicMapOrig = retinotopicMapFromIndividualMaps(mapORIG);
%% ============================= TSCA ===================================
noiseNew.time = eye(T)/T; noiseNew.space = []; mapTSCA = zeros(m,m,length(signals));
noise2New.time = createToeplitz(3,0.1,1,1,T); noise2New.space = [];
noise3New.time = createToeplitz(1,0.1,1,1,T); noise3New.space = [];
for i=1:length(signals)
    sig.time = normpdf(0:0.1:(T-1)/10,stim_time(i),1);
    [projected,components,D,Alpha,output] = tscaFunc(Z - mean(Z(:)),sig,[noiseNew noise2New noise3New],[1 -0.25*ones(1,3)],100,1);
%     tscaAnalyze(output,3,[],0,T);
    [~,I] = max(corr(abs(projected'),sig.time')); % get the index of component with highest correlation to original signal
    mapTSCA(:,:,i) = abs(reshape(components(:,I),m,m));
    [tscaMSE(i),tscaPSNR(i),tscaCNR(i),tscaMSSIM(i),tscaCorr(i),tscaCP(i)] = getPerformance(mapTSCA(:,:,i),signals(i).space,signals(i).space~=DC,signals(i).space==DC);
end
retinotopicMapTSCA = retinotopicMapFromIndividualMaps(mapTSCA);
% [maxmap,maxind] = max(map,[],3); maxmap = maxmap(:); maxind = maxind(:);
% map22 = colors(maxind,:); map22(maxmap(:)< prctile(maxmap(:),70),:) = repmat([0 0 0],length(find(maxmap(:)< prctile(maxmap(:),70))),1);
% retinotopicMap = hsv2rgb(maxind/length(unique(maxind)),ones(size(maxind)),(maxmap-min(maxmap))./(max(maxmap)-min(maxmap)));
% figure;imagesc(reshape(retinotopicMap,40,40,3));
%% ========================== Nadav's method ===========================
for i=1:length(signals)
    thisSig = Z(:,stim_time(i)*10-99:stim_time(i)*10+100);
    [ind_mx1 w_mx i_rx w_rx]=xplmts(thisSig,[],mean(thisSig(:))+0.2*std(thisSig(:)),[0.9 1.1],[mean(thisSig(:))+1*std(thisSig(:)) 30], 10);
    mapN(:,:,i) = reshape(ind_mx1,40,40);
%     figure;imagesc(reshape(mapN(:,:,i),40,40));colormap(gray);
    [nadavMSE(i),nadavPSNR(i),nadavCNR(i),nadavMSSIM(i),nadavCorr(i),nadavCP(i)] = getPerformance(mapN(:,:,i),signals(i).space,signals(i).space~=DC,signals(i).space==DC);
end
retinotopicMapNADAV = retinotopicMapFromIndividualMaps(mapN);
%% ========================== T_max method =============================
refff = normpdf(0:0.1:(T-1)/10,10,1);
for i=1:size(Z,1)
    [rtemp,lags] = xcorr(Z(i,:),refff);
    [r(i),I] = max(rtemp);
    tmax(i) = lags(I);
end
tmax(r<mean(r)+0.2*std(r)) = -100;
A = tmax;
B = [-100, 0:200:800];  
[~,I] = min(abs(bsxfun(@minus,A,B')));
Anew = B(I);
for i=1:length(signals)
    temp = zeros(size(Anew))';
    temp(Anew == B(i+1)) = 1;
    temp(temp==1) = r(temp==1);
    mapT(:,:,i) = reshape(temp,40,40);
    [T_MSE(i),T_PSNR(i),T_CNR(i),T_MSSIM(i),T_corr(i),T_CP(i)] = getPerformance(mapT(:,:,i),signals(i).space,signals(i).space~=DC,signals(i).space==DC);
end
retinotopicMapTmax = retinotopicMapFromIndividualMaps(mapT);
% tmax2 = Anew';
% tmax2(tmax2>=0) = (tmax2(tmax2>=0)./100+2)./2; tmax22 = zeros(length(tmax2),3);% transform 0,200,400,600,800 to 1,2,3,4,5 
% tmax22(tmax2>0,:) = colors(tmax2(tmax2>0),:);
% figure;imagesc(reshape(tmax22,40,40,3));
% close all;
runSummary{1,k} = [origMSE;origPSNR;origCNR;origMSSIM;origCorr;origCP];
runSummary{2,k} = [nadavMSE;nadavPSNR;nadavCNR;nadavMSSIM;nadavCorr;nadavCP];
runSummary{3,k} = [tscaMSE;tscaPSNR;tscaCNR;tscaMSSIM;tscaCorr;tscaCP];
runSummary{4,k} = [T_MSE;T_PSNR;T_CNR;T_MSSIM;T_corr;T_CP];
end
%% Compare results (all runs)
ORIG = cat(3,runSummary{1,:});
NADAV = cat(3,runSummary{2,:}); NADAV(isnan(NADAV)) = 0;
TSCA = cat(3,runSummary{3,:});
Tmax = cat(3,runSummary{4,:}); Tmax(isinf(Tmax)) = 100; Tmax(isnan(Tmax)) = 0;
Title = {'MSE','PSNR','CNR','MSSIM','Corr','CP'};
Title2 = {'TSCA','T max','Nadav','Original'};
retMaps = {retinotopicMapTSCA,retinotopicMapTmax,retinotopicMapNADAV,retinotopicMapOrig};
% figure(kk);suptitle(['Noise STD = ' num2str(noiseSig(kk))]); 
figure(kk*10);suptitle(['Theoretical SNR = ' num2str(10*log10(0.3/noiseSig(kk)))]);
figure(kk*100);suptitle(['Theoretical SNR = ' num2str(10*log10(0.3/noiseSig(kk)))]);
for i=1:6 % iterate the 6 performance measures
%     figure(kk);
%     subplot(2,3,i); 
%     hold on;
%     plot(mean(TSCA(i,:,:),3),'*r','MarkerSize',10);
%     plot(mean(Tmax(i,:,:),3),'*g','MarkerSize',10);
%     plot(mean(NADAV(i,:,:),3),'*k','MarkerSize',10);
%     plot(mean(ORIG(i,:,:),3),'*b','MarkerSize',10);
%     legend('tsca','T','Nadav','orig'); title([Title{i} ' between reconstructed and original image']); ylabel(Title{i}); xlabel('sig #');
    figure(kk*10);
    subplot(2,3,i) 
    boxplot([squeeze(mean(TSCA(i,:,:),2)) squeeze(mean(Tmax(i,:,:),2)) squeeze(mean(NADAV(i,:,:),2)) squeeze(mean(ORIG(i,:,:),2))],char({'tsca';'T';'Nadav';'orig'}));
    title(Title{i});
    if i<=4
       figure(kk*100);
       subplot(2,2,i)
       imagesc(reshape(retMaps{i},40,40,3)); title([Title2{i} ' Retinotopic Map']);
       
    end
end
end