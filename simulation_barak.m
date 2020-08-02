%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% SIMULATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath("C:\Users\orica\Dropbox\fcns_and_decript");
% addpath('C:\Users\orica\Dropbox\master degree\codes');
global fs 
T = 1000; fs = 100;
P = 1600;
m = sqrt(P);
signals = struct;
noises = struct;
stim_time = 10:20:90;
noiseSig = 0.02; % STD of the noise, change to 0.3 for theoretical SNR = 0 
%% construct signals
[I,J] = ndgrid(1:m,1:m); r = 4; locs = [0.25 0.25; 0.25 0.75; 0.5 0.5; 0.75 0.25; 0.75 0.75]; 
DC = 0.02; 
figure;ind2plot = [1,5,9,13,17];
for i=1:5
    signals(i).time = 1*normpdf(0:0.1:(T-1)/10,stim_time(i),1);
    signals(i).space = imgaussfilt(double((I-m*locs(i,1)).^2+(J-m*locs(i,2)).^2<r^2),sqrt(r))+DC; % what is larger than r^2 is 1 (white), what is smaller is 0 (black)
    subplot(3,6,ind2plot(i)) 
    imagesc(signals(i).space); colormap(gray);
    title(['sig. ' num2str(i) ' - spatial']);
    subplot(3,6,ind2plot(i)+1)
    plot(signals(i).time);
    title(['sig. ' num2str(i) ' - temporal']);
    xticks([0 100:200:T]);
end
%% construct noise
t = linspace(0,(T-1)./fs,T);
[I1,I2] = ndgrid([repmat(linspace(0,2*pi,m/2),1,2)]',[repmat(linspace(0,2*pi,m/2),1,2)]);
noises(1).time = normrnd(0,noiseSig,T,1);
noises(1).space = cos(I1);
noises(2).time = noiseSig*cos(2*pi*1.*t')+normrnd(0,noiseSig,T,1);
noises(2).space = normrnd(0,1,m,m);
noises(3).time = noiseSig*cos(2*pi*3.*t')+normrnd(0,noiseSig,T,1);
noises(3).space = cos(I2);
figure; ind2plot = [1,3,5];
for i=1:length(noises)
    subplot(3,2,ind2plot(i))
    imagesc(noises(i).space); colormap(gray);
    title(['noise ' num2str(i) ' - spatial']);
    subplot(3,2,ind2plot(i)+1)
    plot(noises(i).time);
    title(['noise ' num2str(i) ' - temporal']);
end
%% construct Z and show 9 frames
Z = zeros(m*m,T); % preallocate memory
for i = 1:T % can "play" with the noises being added
    Z(:,i) = reshape( ...
        signals(1).time(i)*signals(1).space+signals(2).time(i)*signals(2).space+signals(3).time(i)*signals(3).space+signals(4).time(i)*signals(4).space+signals(5).time(i)*signals(5).space+...
        noises(1).time(i)*noises(1).space+...
        noises(2).time(i)*noises(2).space+... 
        noises(3).time(i)*noises(3).space,...
        [],1);
end
figure; suptitle('the recorded signal (Z) at random frames');
ind2show = sort([100;502;899;randi(1000,7,1)]);
for i=1:length(ind2show)
    subplot(2,5,i)
    imshow(reshape(Z(:,ind2show(i)),m,[])); colormap(gray);
    title(['Frame #: ' num2str(ind2show(i))]);
end
%% ============================= TSCA ===================================
map = zeros(m,m,length(signals));
noiseNew.time = eye(T)/T; noiseNew.space = [];  % autocorrelatoin of white noise
noise2New.time = createToeplitz(3,0.1,1,1,T); noise2New.space = []; % autocorrelation of cos 3[Hz] with noise
noise3New.time = createToeplitz(1,0.1,1,1,T); noise3New.space = []; % autocorrelation of cos 1[Hz] with noise
for i=1:length(signals) % iterate the signals
    sig.time = normpdf(0:0.1:(T-1)/10,stim_time(i),1); % constuct the time course of this signal (tscaFunc calculates the autocorrelation matrix if given a vector by E[vv'])
    [projected,components,D,Alpha,output] = tscaFunc(Z - mean(Z(:)),sig,[noiseNew noise2New noise3New],[1 -0.25*ones(1,3)],100,1); % call function: Raw Input,Signal,Noises,Gammas,# components to project,Perform PCA [0,1]
%     tscaAnalyze(output,3,[],0,T); % this function shows the first n components, the eigenvalues and the autocorrelation matrices 
    [~,I] = max(corr(abs(projected'),sig.time')); % get the index of component with highest correlation to original signal (not always the first one...)
    map(:,:,i) = abs(reshape(components(:,I),m,m));
    [tscaMSE(i),tscaPSNR(i),tscaCNR(i),tscaMSSIM(i),tscaCorr(i),tscaCP(i)] = getPerformance(map(:,:,i),signals(i).space,signals(i).space~=DC,signals(i).space==DC); % get performance measures for each ind. map
end
retinotopicMapTSCA = retinotopicMapFromIndividualMaps(map,1); % create retinotopic map using hsv2rgb
