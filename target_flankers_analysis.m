%% After obtaining maps from VSDI_postAnalysis_main
%% imshow the general responses obtained by tsca&glm
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
imf2(rshp(tempbrn)); hold on; imf2(avgTargetSpatial,0.75);
title('target');
subplot(1,3,2)
avgFlankersSpatial = mean(mapTSCA(:,:,4:6),3);
imf2(rshp(tempbrn)); hold on; imf2(avgFlankersSpatial,0.75);
title('flankers');
subplot(1,3,3)
avgBothSpatial = mean(mapTSCA(:,:,7:9),3);
imf2(rshp(tempbrn)); hold on; imf2(avgBothSpatial,0.75);
title('combined');
%% get target and flankers ROI
[row,col] = find( avgTargetSpatial>=0.75);
targetROI = zeros(sz); 
for i = 1:length(row)
    targetROI(row(i),col(i)) = 1; 
end
figure;imf2(targetROI);
indTargetROI = sub2ind(sz,row,col);

[row,col] = find( avgFlankersSpatial>=0.75);
flankersROI = zeros(sz); 
for i = 1:length(row)
    flankersROI(row(i),col(i)) = 1; 
end
figure;imf2(flankersROI);
indFlankersROI = sub2ind(sz,row,col);
%% calculate average temporal repsones
t = linspace(0,size(Z,2)/100,size(Z,2));
figure;plot(t,mean(ZZ(indTargetROI,:))); 
xlabel('time [sec]'); ylabel('STD'); title('target ROI over time');
avgTarget = mean([ZZ(indTargetROI,1:100);ZZ(indTargetROI,201:300)]);avgTarget = avgTarget - avgTarget(1); 
avgFlankers = mean([ZZ(indTargetROI,301:400);ZZ(indTargetROI,401:500);ZZ(indTargetROI,501:600)]);avgFlankers = (avgFlankers - avgFlankers(1))/1.5; 
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
avgTargetF = mean([ZZ(indFlankersROI,1:100);ZZ(indFlankersROI,201:300)]);avgTargetF = (avgTargetF - avgTargetF(1))/1.5; 
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
plot(t(1:100),avgTarget);
hold on;
plot(t(1:100),avgFlankers);
plot(t(1:100),avgCombined);
legend('target','flankers','combined');xlabel('time [sec]'); ylabel('STD');
title('average temporal response in target ROI');
%% show target ROI target + flankers and combined overlayed
figure;
plot(t(1:100),avgTarget+avgFlankers);
hold on;
plot(t(1:100),avgCombined);
legend('linear summation','combined');xlabel('time [sec]'); ylabel('STD');
title('average temporal response in target ROI');
%% show target ROI target, flankers, linear summation and combined overlayed
figure;
plot(t(1:100),avgTarget);
hold on;
plot(t(1:100),avgFlankers);
plot(t(1:100),avgTarget+avgFlankers);
plot(t(1:100),avgCombined);
legend('target','flankers','linear summation','combined');xlabel('time [sec]'); ylabel('STD');
title('average temporal response in target ROI');