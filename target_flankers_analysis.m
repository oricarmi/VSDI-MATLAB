%% After obtaining maps from VSDI_postAnalysis_main
% load("C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\2021.01.18\mapTSCA.mat");
close all;
%% imshow the general responses obtained by tsca&glm
ZZ3d = rshp(ZZ);
figure;
tempbrn = repmat(reshape((brn-min(brn))./(max(brn)-min(brn)),[],1),1,1,3);
for i=1:9
    subplot(3,3,i)
    imf2(rshp(tempbrn)); hold on; imf2(mapTSCA(:,:,i),0.75);
    if i<4
        title('target');
    elseif i<7
        title('flankers');
    else 
        title('combined');
    end
end
%% imshow average target, flankers, combined
figure;
subplot(1,3,1);
avgTargetSpatial = mean(mapTSCA(:,:,1:3),3);
imf2(rshp(tempbrn)); hold on; imf2(avgTargetSpatial,0.85);
title('target');
subplot(1,3,2)
avgFlankersSpatial = mean(mapTSCA(:,:,4:6),3);
imf2(rshp(tempbrn)); hold on; imf2(avgFlankersSpatial,0.85);
title('flankers');
subplot(1,3,3)
avgCombinedSpatial = mean(mapTSCA(:,:,7:9),3);
imf2(rshp(tempbrn)); hold on; imf2(avgCombinedSpatial,0.85);
title('combined');
figure;
[~,r] = retinotopicMapFromIndividualMaps_target_flanker(cat(3,avgTargetSpatial,avgFlankersSpatial),1,'',97);
%% get target and flankers ROI
ZZ3d = rshp(ZZ);
if input('choose manually?')
    figure;targetROI_manual = roipoly(mean(ZZ3d(:,:,07:15),3)+mean(ZZ3d(:,:,107:115),3)+mean(ZZ3d(:,:,207:215),3));
    [row,col] = find( targetROI_manual);
    indTargetROI = sub2ind(sz,row,col);
    figure;flankersROI1 = roipoly(mean(ZZ3d(:,:,307:315),3)+mean(ZZ3d(:,:,407:415),3)+mean(ZZ3d(:,:,507:515),3));
    figure;flankersROI2 = roipoly(mean(ZZ3d(:,:,307:315),3)+mean(ZZ3d(:,:,407:415),3)+mean(ZZ3d(:,:,507:515),3));
    [row,col] = find( flankersROI1 | flankersROI2);
    flankersROI_manual = flankersROI2 | flankersROI1;
    indFlankersROI = sub2ind(sz,row,col);
    figure;
    subplot 121
    imf2(targetROI_manual);title('target ROI');
    subplot 122
    imf2(flankersROI_manual);title('flankers ROI');
else
    [row,col] = find( imopen(double(avgTargetSpatial>=0.73),strel('disk',2)) );
    targetROI = zeros(sz); 
    for i = 1:length(row)
        targetROI(row(i),col(i)) = 1; 
    end
%     targetROI = imopen(targetROI,strel('disk',5));
    indTargetROI = sub2ind(sz,row,col);
    [row,col] = find( imopen(double(avgFlankersSpatial>=0.83),strel('disk',2)) );
    flankersROI = zeros(sz); 
    for i = 1:length(row)
        flankersROI(row(i),col(i)) = 1; 
    end
    indFlankersROI = sub2ind(sz,row,col);
    figure;
    subplot 121
    imf2(targetROI);title('target ROI');
    subplot 122
    imf2(flankersROI);title('flankers ROI');
end
%% calculate average temporal repsones
t = linspace(0,size(Z,2)/100,size(Z,2));
figure;plot(t,mean(ZZ(indTargetROI,:))); 
xlabel('time [sec]'); ylabel('STD'); title('target ROI over time');
avgTarget = mean([ZZ(indTargetROI,1:100);ZZ(indTargetROI,101:200);ZZ(indTargetROI,201:300)]);avgTarget = (avgTarget - avgTarget(1))*1.2; 
avgFlankers = mean([ZZ(indTargetROI,301:400);ZZ(indTargetROI,401:500);ZZ(indTargetROI,501:600)]);avgFlankers = (avgFlankers - avgFlankers(1))/2; 
avgCombined = mean([ZZ(indTargetROI,601:700);ZZ(indTargetROI,701:800);ZZ(indTargetROI,801:900)]);avgCombined = avgCombined - avgCombined(1); 
%% show temporal evolution of target ROI
figure;
suptitle('target ROI');
subplot(2,3,1);
plot(t(1:300),mean(ZZ(indTargetROI,1:300)));
title('target');xlabel('time [sec]'); ylabel('STD');
subplot(2,3,4);
% plot(mean([Z(indROI,1:100);Z(indROI,101:200);Z(indROI,201:300)]))
plot(t(1:100),avgTarget)
title('target');xlabel('time [sec]'); ylabel('STD');

subplot(2,3,2)
plot(t(1:300),mean(ZZ(indTargetROI,301:600)));
title('flankers');xlabel('time [sec]'); ylabel('STD');
subplot(2,3,5);
plot(t(1:100),avgFlankers)
title('flankers');xlabel('time [sec]'); ylabel('STD');

subplot(2,3,3)
plot(t(1:300),mean(ZZ(indTargetROI,601:900)));
title('combined');xlabel('time [sec]'); ylabel('STD');
subplot(2,3,6);
plot(t(1:100),avgCombined)
title('combined');xlabel('time [sec]'); ylabel('STD');

%% show temporal evolution of flankers ROI
t = linspace(0,size(Z,2)/100,size(Z,2));
figure;plot(t,mean(ZZ(indFlankersROI,:)));
xlabel('time [sec]'); ylabel('STD'); title('flankers ROI over time');
avgTargetF = mean([ZZ(indFlankersROI,1:100);ZZ(indFlankersROI,201:300)]);avgTargetF = (avgTargetF - avgTargetF(1)); 
avgFlankersF = mean([ZZ(indFlankersROI,301:400);ZZ(indFlankersROI,401:500);ZZ(indFlankersROI,501:600)]);avgFlankersF = avgFlankersF - avgFlankersF(1); 
avgCombinedF = mean([ZZ(indFlankersROI,601:700);ZZ(indFlankersROI,701:800);ZZ(indFlankersROI,801:900)]);avgCombinedF = avgCombinedF - avgCombinedF(1); 

figure;
suptitle('flankers ROI');
subplot(2,3,1);
plot(t(1:300),mean(ZZ(indFlankersROI,1:300)));
title('target');xlabel('time [sec]'); ylabel('STD');
subplot(2,3,4);
% plot(mean([Z(indROI,1:100);Z(indROI,101:200);Z(indROI,201:300)]))
plot(t(1:100),avgTargetF)
title('target');xlabel('time [sec]'); ylabel('STD');

subplot(2,3,2)
plot(t(1:300),mean(ZZ(indFlankersROI,301:600)));
title('flankers');xlabel('time [sec]'); ylabel('STD');
subplot(2,3,5);
plot(t(1:100),avgFlankersF)
title('flankers');xlabel('time [sec]'); ylabel('STD');

subplot(2,3,3)
plot(t(1:300),mean(ZZ(indFlankersROI,601:900)));
title('combined');xlabel('time [sec]'); ylabel('STD');
subplot(2,3,6);
plot(t(1:100),avgCombinedF)
title('combined');xlabel('time [sec]'); ylabel('STD');
%% show target ROI target, flankers and combined overlayed
figure;
plot(t(1:100),avgFlankers);
xlabel('time [sec]'); ylabel('STD');title('average temporal response in target ROI');
xline(0.07,'--r'); text(0.072,-0.5,'70 [msec]','color','r');xline(0.35,'--r'); text(0.352,-0.5,'350 [msec]','color','r');
% hold on;
% plot(t(1:100),avgFlankers);
% plot(t(1:100),avgCombined);
% legend('target','flankers','combined');
% title('average temporal response in target ROI');
% %% show target ROI target + flankers and combined overlayed
% figure;
% plot(t(1:100),avgTarget+avgFlankers);
% hold on;
% plot(t(1:100),avgCombined);
% legend('linear summation','combined');xlabel('time [sec]'); ylabel('STD');
% title('average temporal response in target ROI');
%% show target ROI target, flankers, linear summation and combined overlayed
figure;
plot(t(1:100),avgTarget);
hold on;
plot(t(1:100),avgFlankers);
plot(t(1:100),avgTarget+avgFlankers);
plot(t(1:100),avgCombined);
legend('target','flankers','linear summation','combined');xlabel('time [sec]'); ylabel('STD');
title('average temporal response in target ROI');
%% Cross correlation between target & flanker and target & combined
% [rTF,lags] = xcorr(avgTarget,avgFlankers,20,'coeff');
% [rTC] = xcorr(avgTarget,avgCombined,20,'coeff');
% [rFC] = xcorr(avgFlankers,avgCombined,20,'coeff');
% figure;
% subplot 311
% plot(lags./fs*1000,rTF);title('cross correlation between target and flankers'); xlabel('lags [msec]'); ylabel('\rho');
% subplot 312
% plot(lags./fs*1000,rTC);title('cross correlation between target and combined'); xlabel('lags [msec]'); ylabel('\rho');
% subplot 313
% plot(lags./fs*1000,rFC);title('cross correlation between flankers and combined'); xlabel('lags [msec]'); ylabel('\rho');
%% upsample signals in target ROI
avgTarget2 = interp(avgTarget,5);
avgFlankers2 = interp(avgFlankers,5);
avgCombined2 = interp(avgCombined,5);
t2 = linspace(0,size(Z,2)/100,size(Z,2)*5);
figure;
plot(t2(1:500),avgTarget2);
hold on;
plot(t2(1:500),avgFlankers2);
plot(t2(1:500),avgTarget2+avgFlankers2);
plot(t2(1:500),avgCombined2);
legend('target','flankers','linear summation','combined');xlabel('time [sec]'); ylabel('STD');
title('avg. temporal res. (upsampled) - target ROI ');
%% Cross correlation between target & flanker and target & combined (upsampled)
[rTF2,lags2] = xcorr(avgTarget2,avgFlankers2,20,'coeff');
[rTC2] = xcorr(avgTarget2,avgCombined2,20,'coeff');
[rFC2] = xcorr(avgFlankers2,avgCombined2,20,'coeff');
rTF2 = fliplr(rTF2);
lagsTime = lags2./(fs*5)*1000;
figure;
subplot 311
plot(lagsTime,rTF2);title('cross correlation between target and flankers'); xlabel('lags [msec]'); ylabel('\rho');%ylim([0.8 0.95]);
hold on;plot(lagsTime(rTF2==max(rTF2)),max(rTF2),'.r','MarkerSize',20);
text(lagsTime(rTF2==max(rTF2)),max(rTF2)-0.05,['max corr. at ' num2str(lagsTime(rTF2==max(rTF2))) ' [msec]'],'FontSize',10);
subplot 312
plot(lagsTime,rTC2);title('cross correlation between target and combined'); xlabel('lags [msec]'); ylabel('\rho');%ylim([0.8 0.95]);
hold on;plot(lagsTime(rTC2==max(rTC2)),max(rTC2),'.r','MarkerSize',20);
text(lagsTime(rTC2==max(rTC2)),max(rTC2)-0.05,['max corr. at ' num2str(lagsTime(rTC2==max(rTC2))) ' [msec]'],'FontSize',10);
subplot 313
plot(lagsTime,rFC2);title('cross correlation between flankers and combined'); xlabel('lags [msec]'); ylabel('\rho');%ylim([0.8 0.95]);
hold on;plot(lagsTime(rFC2==max(rFC2)),max(rFC2),'.r','MarkerSize',20);
text(lagsTime(rFC2==max(rFC2)),max(rFC2)-0.01,['max corr. at ' num2str(lagsTime(rFC2==max(rFC2))) ' [msec]'],'FontSize',10);
%% upsample signals in flanker ROI
avgTargetF2 = interp(avgTargetF,5);
avgFlankersF2 = interp(avgFlankersF,5);
avgCombinedF2 = interp(avgCombinedF,5);
t2 = linspace(0,size(Z,2)/100,size(Z,2)*5);
figure;
plot(t2(1:500),avgTargetF2);
hold on;
plot(t2(1:500),avgFlankersF2);
plot(t2(1:500),avgTargetF2+avgFlankersF2);
plot(t2(1:500),avgCombinedF2);
legend('target','flankers','linear summation','combined');xlabel('time [sec]'); ylabel('STD');
title('average temporal response in flankers ROI');

%% Cross correlation between target & flanker and target & combined (upsampled) - Flankers ROI
[rTF22,lags22] = xcorr(avgTargetF2,avgFlankersF2,20,'coeff');
[rTC22] = xcorr(avgTargetF2,avgCombinedF2,20,'coeff');
[rFC22] = xcorr(avgFlankersF2,avgCombinedF2,20,'coeff');
lagsTime22 = lags22./(fs*5)*1000;
figure;
subplot 311
plot(lagsTime22,rTF22);title('cross correlation between target and flankers'); xlabel('lags [msec]'); ylabel('\rho');%ylim([0.8 0.95]);
hold on;plot(lagsTime22(rTF22==max(rTF22)),max(rTF22),'.r','MarkerSize',20);
text(lagsTime22(rTF22==max(rTF22)),max(rTF22)-0.05,['max corr. at ' num2str(lagsTime22(rTF22==max(rTF22))) ' [msec]'],'FontSize',10);
subplot 312
plot(lagsTime22,rTC22);title('cross correlation between target and combined'); xlabel('lags [msec]'); ylabel('\rho');%ylim([0.8 0.95]);
hold on;plot(lagsTime22(rTC22==max(rTC22)),max(rTC22),'.r','MarkerSize',20);
text(lagsTime22(rTC22==max(rTC22)),max(rTC22)-0.05,['max corr. at ' num2str(lagsTime22(rTC22==max(rTC22))) ' [msec]'],'FontSize',10);
subplot 313
plot(lagsTime22,rFC22);title('cross correlation between flankers and combined'); xlabel('lags [msec]'); ylabel('\rho');%ylim([0.8 0.95]);
hold on;plot(lagsTime22(rFC22==max(rFC22)),max(rFC22),'.r','MarkerSize',20);
text(lagsTime22(rFC22==max(rFC22)),max(rFC22)-0.02,['max corr. at ' num2str(lagsTime22(rFC22==max(rFC22))) ' [msec]'],'FontSize',10);
%% derivative maps - all brain
ZZ3d = reshape(ZZ,[],100,9);
derZZ = zeros(size(ZZ3d));
for i=1:size(ZZ3d,2)
    if i<100
        derZZ(:,i,:) = ZZ3d(:,i+1,:)-ZZ3d(:,i,:);
    else
        derZZ(:,i,:) = ZZ3d(:,i,:);
    end
end
avgDerZZ_target = mean(derZZ(:,:,1:3),3);
avgDerZZ_flankers = mean(derZZ(:,:,4:6),3);
avgDerZZ_combined = mean(derZZ(:,:,7:8),3);
% figure(101);suptitle('target derivative maps');
% figure(102);suptitle('flankers derivative maps');
% figure(103);suptitle('combined derivative maps');
% time_after_onset = 10:20:500; % in msec
% for i=1:length(time_after_onset)
%     thisTime = time_after_onset(i);
%     figure(101);
%     subplot(sqrt(length(time_after_onset)),sqrt(length(time_after_onset)),i);
%     imf2(rshp(avgDerZZ_target(:,thisTime/10)));title([ num2str(thisTime) 'msec after stimulus onset']);
%     figure(102);
%     subplot(sqrt(length(time_after_onset)),sqrt(length(time_after_onset)),i);
%     imf2(rshp(avgDerZZ_flankers(:,thisTime/10)));title([ num2str(thisTime) 'msec after stimulus onset']);    
%     figure(103);
%     subplot(sqrt(length(time_after_onset)),sqrt(length(time_after_onset)),i);
%     imf2(rshp(avgDerZZ_combined(:,thisTime/10)));title([ num2str(thisTime) 'msec after stimulus onset']);
% end
%% derivative maps - target ROI
time_after_onset = 10:20:500; % in msec
ZZ3d = reshape(ZZ,[],100,9);
derZZ = zeros(size(ZZ3d));
filtWidth = 3;
filtSigma = 1;
imageFilter=fspecial('gaussian',filtWidth,filtSigma);
not_indTargetROI = setdiff(1:size(ZZ,1),indTargetROI);
for i=1:size(ZZ3d,2)
%     if i<100
%         derZZ(indTargetROI,i,:) = ZZ3d(indTargetROI,i+1,:)-ZZ3d(indTargetROI,i,:);
%     else
%         derZZ(indTargetROI,i,:) = ZZ3d(indTargetROI,i,:);
%     end
    if i>1 && i<100
        derZZ(indTargetROI,i,:) = ZZ3d(indTargetROI,i+1,:)-ZZ3d(indTargetROI,i-1,:);
    else
        derZZ(indTargetROI,i,:) = ZZ3d(indTargetROI,i,:);
    end
end
avgDerZZ_target = reshape(imgaussfilt(rshp(mean(derZZ(:,:,1:3),3)),0.6),[],100); 
avgDerZZ_target(not_indTargetROI,:) = nan; 
% avgDerZZ_target = (max(avgDerZZ_target(:))-min(avgDerZZ_target(:)))*MinMaxNorm(avgDerZZ_target) + min(avgDerZZ_target(:));
avgDerZZ_target1 = 2*( (avgDerZZ_target - repmat(min(avgDerZZ_target,[],1),size(avgDerZZ_target,1),1))./(repmat(max(avgDerZZ_target,[],1),size(avgDerZZ_target,1),1)-repmat(min(avgDerZZ_target,[],1),size(avgDerZZ_target,1),1)) );
figure(201);suptitle('target derivative maps - all same scale');
figure(2011);suptitle('target derivative maps - normalized (same scale)');
avgDerZZ_flankers = reshape(imgaussfilt(rshp(mean(derZZ(:,:,4:6),3)),0.4),[],100);
avgDerZZ_flankers(not_indTargetROI,:) = nan; 
% avgDerZZ_flankers = 2*MinMaxNorm(avgDerZZ_flankers) + min(avgDerZZ_flankers(:));
avgDerZZ_flankers1 = 2*( (avgDerZZ_flankers - repmat(min(avgDerZZ_flankers,[],1),size(avgDerZZ_flankers,1),1))./(repmat(max(avgDerZZ_flankers,[],1),size(avgDerZZ_flankers,1),1)-repmat(min(avgDerZZ_flankers,[],1),size(avgDerZZ_flankers,1),1)) );
avgDerZZ_combined = reshape(imgaussfilt(rshp(mean(derZZ(:,:,7:9),3)),0.4),[],100); 
avgDerZZ_combined(not_indTargetROI,:) = nan; 
avgDerZZ_combined = 255*MinMaxNorm(avgDerZZ_combined);
avgDerZZ_combined1 = 255*( (avgDerZZ_combined - repmat(min(avgDerZZ_combined,[],1),size(avgDerZZ_combined,1),1))./(repmat(max(avgDerZZ_combined,[],1),size(avgDerZZ_combined,1),1)-repmat(min(avgDerZZ_combined,[],1),size(avgDerZZ_combined,1),1)) );
figure(202);suptitle('flankers derivative maps - all same scale');
figure(2021);suptitle('flankers derivative maps - normalized (same scale)');
figure(203);suptitle('combined derivative maps - all same scale');
figure(2031);suptitle('combined derivative maps - normalized (same scale)');
% cAxis = [prctile([avgDerZZ_target(:);avgDerZZ_flankers(:);avgDerZZ_combined(:)],1),prctile([avgDerZZ_target(:);avgDerZZ_flankers(:);avgDerZZ_combined(:)],98)];
cAxis1 = [prctile([avgDerZZ_target(:)],1),prctile([avgDerZZ_target(:)],99)];
cAxis2 = [prctile([avgDerZZ_flankers(:)],1),prctile([avgDerZZ_flankers(:)],99)];
xLim = [1.55,3.45];yLim = [2,4.8];
for i=1:length(time_after_onset)
    thisTime = time_after_onset(i);
    figure(201);
    subplot(sqrt(length(time_after_onset)),sqrt(length(time_after_onset)),i);
    imf2(rshp(avgDerZZ_target(:,thisTime/10)));title([ num2str(thisTime) 'msec']);
    xlabel('');set(gca,'xticklabel',[]);set(gca,'xtick',[]);ylabel('');set(gca,'yticklabel',[]);set(gca,'ytick',[]);
    xlim(xLim);ylim(yLim);
    caxis(cAxis1);
    if i==1
        colorbar;
        rectangle('Position',[xLim(1)+0.1,yLim(2)-0.2,1,0.1],'FaceColor',[1 1 1],'EdgeColor',[1 1 1]);
        text(xLim(1)+0.3,yLim(2)-0.25,'1mm','color',[1 1 1],'fontsize',8);
    end
    figure(2011);
    subplot(sqrt(length(time_after_onset)),sqrt(length(time_after_onset)),i);
    imf2(rshp(avgDerZZ_target1(:,thisTime/10)));title([ num2str(thisTime) 'msec']);
    xlabel('');set(gca,'xticklabel',[]);set(gca,'xtick',[]);ylabel('');set(gca,'yticklabel',[]);set(gca,'ytick',[]);
    xlim(xLim);ylim(yLim);
    if i==1
        colorbar;
        rectangle('Position',[xLim(1)+0.1,yLim(2)-0.2,1,0.1],'FaceColor',[1 1 1],'EdgeColor',[1 1 1]);
        text(xLim(1)+0.3,yLim(2)-0.25,'1mm','color',[1 1 1],'fontsize',8);
    end
    figure(202);
    subplot(sqrt(length(time_after_onset)),sqrt(length(time_after_onset)),i);
    imf2(rshp(avgDerZZ_flankers(:,thisTime/10)));title([ num2str(thisTime) 'msec after stimulus onset']);  
    xlabel('');set(gca,'xticklabel',[]);set(gca,'xtick',[]);ylabel('');set(gca,'yticklabel',[]);set(gca,'ytick',[]);
    xlim(xLim);ylim(yLim);
    caxis(cAxis2);
    if i==1
        colorbar;
        rectangle('Position',[2.9,4.8,1,0.1],'FaceColor',[1 1 1],'EdgeColor',[1 1 1]);
        text(2.9,4.65,'1mm','color',[1 1 1],'fontsize',8);
    end
    figure(2021);
    subplot(sqrt(length(time_after_onset)),sqrt(length(time_after_onset)),i);
    imf2(rshp(avgDerZZ_flankers1(:,thisTime/10)));title([ num2str(thisTime) 'msec']);
    xlabel('');set(gca,'xticklabel',[]);set(gca,'xtick',[]);ylabel('');set(gca,'yticklabel',[]);set(gca,'ytick',[]);
    xlim(xLim);ylim(yLim);
    if i==1
        colorbar;
        rectangle('Position',[2.8,4.8,1,0.1],'FaceColor',[1 1 1],'EdgeColor',[1 1 1]);
        text(2.9,4.65,'1mm','color',[1 1 1],'fontsize',8);
    end
%     figure(203);
%     subplot(sqrt(length(time_after_onset)),sqrt(length(time_after_onset)),i);
%     imf2(rshp(avgDerZZ_combined(:,thisTime/10)));title([ num2str(thisTime) 'msec after stimulus onset']);
%     xlabel('');set(gca,'xticklabel',[]);set(gca,'xtick',[]);ylabel('');set(gca,'yticklabel',[]);set(gca,'ytick',[]);
%     xlim(xLim);ylim(yLim);
%     caxis(cAxis);
%     if i==1
%         colorbar;
%         rectangle('Position',[2.9,4.8,1,0.1],'FaceColor',[1 1 1],'EdgeColor',[1 1 1]);
%         text(2.9,4.65,'1mm','color',[1 1 1],'fontsize',8);
%     end
%     figure(2031);
%     subplot(sqrt(length(time_after_onset)),sqrt(length(time_after_onset)),i);
%     imf2(rshp(avgDerZZ_combined1(:,thisTime/10)));title([ num2str(thisTime) 'msec']);
%     xlabel('');set(gca,'xticklabel',[]);set(gca,'xtick',[]);ylabel('');set(gca,'yticklabel',[]);set(gca,'ytick',[]);
%     xlim(xLim);ylim(yLim);
%     if i==1
%         colorbar;
%         rectangle('Position',[2.8,4.8,1,0.1],'FaceColor',[1 1 1],'EdgeColor',[1 1 1]);
%         text(2.9,4.65,'1mm','color',[1 1 1],'fontsize',8);
%     end
end

%% derivative maps - flankers ROI
ZZ3d = reshape(ZZ,[],100,9);
derZZ = zeros(size(ZZ3d));
for i=1:size(ZZ3d,2)
    if i<100
        derZZ(indFlankersROI,i,:) = ZZ3d(indFlankersROI,i+1,:)-ZZ3d(indFlankersROI,i,:);
    else
        derZZ(indFlankersROI,i,:) = ZZ3d(indFlankersROI,i,:);
    end
end
avgDerZZ_target = mean(derZZ(:,:,1:3),3);avgDerZZ_target(~avgDerZZ_target) = nan;
avgDerZZ_flankers = mean(derZZ(:,:,4:6),3);avgDerZZ_flankers(~avgDerZZ_flankers) = nan;
avgDerZZ_combined = mean(derZZ(:,:,7:8),3);avgDerZZ_combined(~avgDerZZ_combined) = nan;
figure(301);suptitle('target derivative maps');
figure(302);suptitle('flankers derivative maps');
figure(303);suptitle('combined derivative maps');
for i=1:length(time_after_onset)
    thisTime = time_after_onset(i);
    figure(301);
    subplot(sqrt(length(time_after_onset)),sqrt(length(time_after_onset)),i);
    imf2(rshp(tempbrn)); hold on;
    imf2(rshp(avgDerZZ_target(:,thisTime/10)));title([ num2str(thisTime) 'msec after stimulus onset']);
    figure(302);
    subplot(sqrt(length(time_after_onset)),sqrt(length(time_after_onset)),i);
    imf2(rshp(tempbrn)); hold on;
    imf2(rshp(avgDerZZ_flankers(:,thisTime/10)));title([ num2str(thisTime) 'msec after stimulus onset']);    
    figure(303);
    subplot(sqrt(length(time_after_onset)),sqrt(length(time_after_onset)),i);
    imf2(rshp(tempbrn)); hold on;
    imf2(rshp(avgDerZZ_combined(:,thisTime/10)));title([ num2str(thisTime) 'msec after stimulus onset']);
end
