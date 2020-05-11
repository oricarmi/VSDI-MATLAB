%% plot all pixels in time heat map
nm={'C', 'DL', 'UR', 'D', 'U', 'DR', 'UL', 'R', 'L'};
figure;
for i=2:10
    y = cfn{i};
    subplot(3,3,i-1);
    imagesc(y,[prctile(y(:),1) prctile(y(:),99)]);
    %c = colorbar('southoutside'); 
    colormap(gca,jet);
    title(nm{i-1});
    xlabel('time [csec]'); ylabel('pixel');
end
%% LMS
Raw_Signal = cfn{2}(44000,:); % center (and baseline in it)
Noise_ref = cfn{1}(44000,:); % baseline of the same pixel
Options = struct('mu',0.1);
Clean_Signal = LMS( Raw_Signal , Noise_ref, Options );
figure;
subplot 211
plot(Raw_Signal);hold on; plot(Noise_ref);
legend('Center original','Baseline');
subplot 212
plot(Clean_Signal);hold on; plot(Noise_ref);
legend('Center cleaned','Baseline');
%%
x = repmat([1:size(brn,2)]',size(brn,1),1);
y = repmat([1:size(brn,1)]',size(brn,2),1);
z = reshape(img{5},[],1);
M = [ x y z];
GMModel = fitgmdist(M,1);
[~,score] = pca(meas,'NumComponents',2);
fittedImage = mvnrnd(GMModel.mu(1:end-1),GMModel.Sigma(1:end-1,1:end-1),numel(z));
%% draw contours
fitted2DGContours = cell(size(fitted2DG));
zfit2 = zeros(size(zfit));
percents = [95.5 96; 96.5 97; 97.5 98; 98.5 99; 99.5 100];
for j= 1:length(fitted2DGContours)
    for i=1:length(percents)
        zfit2(zfit<prctile(reshape(zfit,[],1),percents(i,2)) & zfit>prctile(reshape(zfit,[],1),percents(i,1))) = 1;
    end
    fitted2DGContours{j} = zfit2;
end
figure; hAxes = gca;
colors = [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1];
MAP = colors(I,:);
colormap( hAxes , [1 1 1; 1 0 0; 0 1 0] )
%%
t = linspace(0,1,100); T = length(t);
figure(4);plot(t,mean(cfn{1}),'linewidth',3);
% [pxx,f] = periodogram(mean(cfn{3}),[],[],100);
% figure;plot(f,pxx);
title('waveform for correlation matrix - noise');xlabel('time [sec]'); ylabel('Amp'); 
figure(5);imagesc([mean(cfn{1})]'*mean(cfn{1})./T); colormap(jet); title('correlation matrix - noise');
figure(6);imagesc(imgaussfilt(reshape(mean(cfn{1},2),size(brn)),3.5)); colormap(jet); title('mean of all frames');
figure(1); suptitle('waveforms for correlation matrices - conditions (signal)');
figure(2); suptitle('correlation matrices - conditions (signal)');
figure(3); suptitle('mean of all frames - conditions (signal)');
index2plot = [11,3,15,7,5,13,1,9]; % 8 locs
for i=3:length(cfn)
    figure(1);
    subplot(3,5,index2plot(i-2));
    plot(t,mean(cfn{i}),'linewidth',2);
    title(lgn(i,:));
    xlabel('time [sec]'); ylabel('Amp'); 
    figure(2);
    subplot(3,5,index2plot(i-2));
    imagesc([mean(cfn{i})]'*mean(cfn{i})./T); colormap(jet);
    title(lgn(i,:));
    figure(3);
    subplot(3,5,index2plot(i-2));
    imagesc(imgaussfilt(reshape(mean(cfn{i},2),size(brn)),3.5)); colormap(jet);
    title(lgn(i,:));
end
%%
Z = cfn{4};
Z = [];
for i=1:length(cfn)
%     if i==2
%         continue;
%     end
    Z = [Z cfn{i}];
end
Z = Z - mean(Z,2);
[U,S,V] = svd(Z,'econ');
[ppx,f] = periodogram(U(1,:),[],[],fs);
[ppx,f] = periodogram(mean(Z),[],[],fs);
[ppx,f] = pmtm(mean(Z),2,pow2(nextpow2(length(mean(Z)))),fs);
figure;
subplot(1,2,1);plot(mean(Z)); xlabel('frame'); ylabel('amp'); title('time course of average signal');
subplot(1,2,2); plot(f,ppx); xlabel('frequency [Hz]');ylabel('amp');title('spectrum of average signal');
figure;
for i=1:49
    subplot(7,7,i)
    imagesc(reshape(U(:,i),size(brn,1),[]));
end
%%
index2plot = [5,7,3,8,2,9,1,6,4];
figure(1000);
for i=1:length(img)
    figure(1000);
    subplot(3,3,index2plot(i));
    imagesc(img{i}); title(lgn(i+1,:)); colormap(jet);caxis([prctile(img{i}(:),50) prctile(img{i}(:),99)]);
end
%%
N = 1024;
t = (0:N-1)';
f = (0:N-1)'/N;
fc = 10/N;
H = @(x) atan(x); % the system. TI, apx linear for small inputs.
a = logspace(-1,2,4);
x = sin(2*pi*fc*t)*a;
X = fft(x);
y = H(x);
Y = fft(y);

for i=1:length(a)
    figure;
    subplot(3,1,1)
    plot(t, y(:,i));
    subplot(3,1,2);
    plot(f, abs(Y(:,i)));
    subplot(3,1,3);
    plot(x(:,i), y(:,i), '.');
end
%%
N = 500;
fs = 100;
t = linspace(0,(N-1)/fs,N)';
x = sin(2*pi*3*t);%+cos(2*pi*3*t);
f = linspace(0,fs,N);X = fft(x);
figure;plot(t,x);figure;plot(f(1:N/2),abs(X(1:N/2)));
y = H(x);
Y = fft(y);

for i=1:length(a)
    figure;
    subplot(3,1,1)
    plot(t, y(:,i));
    subplot(3,1,2);
    plot(f(1:N/2), abs(Y(1:N/2,i)));
    subplot(3,1,3);
    plot(x(:,i), y(:,i), '.');
end
%%
figure;
for i=110:145
    subplot(6,6,i-109);
    imagesc(imgaussfilt(reshape(Z(:,i),size(brn0,1),[]),2.5));colormap(jet);
    title(num2str(i));
end
%%
img2{6}(:,200:end) = 0;
img2{8}(:,200:end) = 0;
img2{7}(116,123) = 1;
figure(500);figure(600);
index2plot = [5,7,3,8,2,9,1,6,4];
for i=1:length(img2)
    disp(i);
    if i~=6 && i~=8 && i~=5
        img2fit = img2{i}; img2fit(img2fit<prctile(img2fit(:),65)) = min(img2fit(:));
    else
        img2fit = img2{i}; img2fit(img2fit<prctile(img2fit(:),85)) = min(img2fit(:));
    end
    [~,fitted2DG{i}] = fmgaussfit(1:size(brn,2),1:size(brn,1),img2fit);
    figure(600);
    subplot(3,3,index2plot(i));
    imagesc(fitted2DG{i}); title(lgn(i+1,:)); colormap(jet);
    figure(500);
    subplot(3,3,index2plot(i));
    imagesc(img2{i}); title(lgn(i+1,:)); colormap(jet); caxis([prctile(img2{i}(:),10) prctile(img2{i}(:),99)]);

end
plotContours9(fitted2DG);
%%
Z = []; % concatenate
for i=1:length(cfn)
%     if i==2 
%         continue
%     end
    Z = [Z cfn{i}];
end
T = size(Z,2);
ZZ = preProcess(Z,T);
[r,c] = ind2sub([size(brn,1),size(brn,2)],27000);
figure;plot(squeeze(mean(mean(ZZ,2),1)));
implay(ZZ,20);
%%
[~,ind] = min(Z,[],2);
figure;imagesc(reshape(min(Z,[],2),size(brn,1),[]));
vec = repmat(170:180,9,1) + repmat([0:50:400]',1,11);
figure;histogram(ZZ(:,:,vec));

%%
index2plot = [5,7,3,8,2,9,1,6,4];
figure(90);
img2 = cell(8,1); fitted2DG = cell(8,1);
for i=2:length(cfn)
    subplot(3,3,index2plot(i-1));
    tempIMG = imgaussfilt(mean(ZZ(:,:,(i-1)*100+6:(i-1)*100+10),3),3);
%     img2{i-1} = tempIMG;
%     [~,fitted2DG{i-1}] = fmgaussfit(1:size(brn,2),1:size(brn,1),tempIMG);
    imagesc(tempIMG); colormap(jet);
    caxis([prctile(tempIMG(:),10) prctile(tempIMG(:),100)]);
end
plotContours9(fitted2DG);
%% constructing retinotopic maps with tmax value
Z = []; % concatenate
for i=4:length(cfn)
%     if i==2 
%         continue
%     end
    Z = [Z cfn{i}];
end
Z = preProcess(Z,size(Z,2));
ZZ = reshape(Z,[],450);
refff = normpdf(0:0.1:(size(Z,2)-1)/10,2,0.5);
for i=1:size(ZZ,1)
    [rtemp,lags] = xcorr(ZZ(i,:),refff);
    [r(i),I] = max(rtemp);
    tmax(i) = lags(I);
end
brnn = squeeze(repmat(reshape((brn-min(brn))./max(brn),[],1),1,1,3));
colors = [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1];
tmaxx = tmax;
tmaxx(r<mean(r)+std(r)) = -100;
A = tmaxx;
B = [-100, 0:50:400];  
[~,I] = min(abs(bsxfun(@minus,A,B')));
Anew = B(I);
tmax2 = reshape(repmat(Anew,1,1,3),[],3);
tmax2(any(mod(tmax2,100)==0,2) & tmax2(:,1)>0,:) = colors(sum(tmax2(any(mod(tmax2,100)==0,2) & tmax2(:,1)>0,:)/100+1,2)/3,:);
tmax2(tmax2(:,1)<0,:) = brnn(tmax2(:,1)<0,:);
for i=1:9
    temp = zeros(size(Anew))';
    temp(Anew == B(i+1)) = 1;
    mapT(:,:,i) = reshape(temp,270,[]);
%     T_MSE(i) = sum((reshape(mapT(:,:,i),[],1) - reshape(signals(i).space,[],1)).^2);
%     T_PSNR(i) = psnr(mapT(:,:,i),signals(i).space);
%     T_corr(i) = corr2(mapT(:,:,i),signals(i).space);
end

figure;imagesc(ump.*[0:size(brn,2)-1]./1000, ump.*[0:size(brn,1)-1]./1000,reshape(tmax2,270,[],3));colormap(jet); 
figure;imagesc(reshape(colors,8,1,3));
%%
for j=1:7
    signalNoise(j).time = [zeros(1,(j)*100) signal2_smooth zeros(1,(length(cfn)-(j+3))*100)];
end
[projected,components,~,~,output] = tscaFunc(Z,signal,[noise1 noise2 signalNoise],[1 -0.01*ones(1,3) -0.01*ones(1,7)] ,800,reduceComp);
tscaAnalyze(output,3,[],1,800);
%%
for i =1:8
    mapTSCA2{i} = MinMaxNorm(mapTSCA{i});
end
retinotopicMapTSCA = retinotopicMapFromIndividualMaps(cat(3,mapTSCA2{:}),2);
%%
Signal = zeros(size(Z));
for j= 1:length(stim_time) % iterate stimulations
    for i=1:size(Z,1) % iterate all pixels
        Signal(i,(stim_time(j)*10-19):(stim_time(j)*10+80)) = GLM_VSDI(Z(i,(stim_time(j)*10-19):(stim_time(j)*10+80)),3,basis(:,1));
    end
end
%%
[Signal1,Signal2,Signal3] = GLM_VSDI(Z,[1 3.3 6.6],theoreticalSig');
%%
theoreticalSig = zeros(8,800);
for i=1:8
    theoreticalSig(i,:) = [zeros(1,(i-1)*100) [basis(:,1)+basis(:,2)]' zeros(1,(8-(i))*100)];
end
maps1 = Z*theoreticalSig';
normalizedZ = (Z-mean(Z,2))./std(Z,[],2);
maps2 = normalizedZ*theoreticalSig';
%%
figure;
index2plot = [11,3,15,7,5,13,1,9]; % 8 locs
cAxis = [prctile(maps2(:),10) prctile(maps2(:),99)];
for i=1:8
   subplot(3,5,index2plot(i));
   imagesc(mapTSCA1(:,:,i));title(lgn(i+2,:)); colormap(jet); caxis([]);
end
for i=1:8
    temp = maps2(:,i);
    map222(:,:,i) = reshape((temp-min(temp))./(max(temp)-min(temp)),size(brn,1),[]);
end
%%
[~,ind1] = sort(squeeze(maps11(116,159,:))')
[~,ind2] = sort(squeeze(maps22(116,159,:))')
[~,indTSCA] = sort(squeeze(mapTSCA1(116,159,:))')
[~,indT] = sort(squeeze(mapT(116,159,:))')
[~,indAOF] = sort(squeeze(mapAOF1(116,159,:))')
%%
plotMaps(mapTSCA,8,'TSCA');
plotMaps(mapAOF,8,'AOF');
plotMaps(mapCorr,8,'Corr');
plotMaps(mapTmax,8,'Tmax');
plotMaps(mapNadav,8,'Nadav');