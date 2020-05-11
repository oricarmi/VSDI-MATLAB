%%
global fs
T = 1000; fs= 100;
P = 1600;
m = sqrt(P);
signal1 = struct;
signal2 = struct;
signal3 = struct;
signal4 = struct;
noise1 = struct;
noise2 = struct;
noise3 = struct;
%% signal temporal and spatial components
I = [repmat(linspace(0,10,T/10),1,10)]';
% signal.time = repelem(abs(normrnd(1,0.0625,10,1)),100).*exp(-1.*I);
% signal.time = ones(T,1);
signal1.time = [zeros(1,20) [exp((-1/20).*[0:(T-21)])]]';
signal2.time = [zeros(1,200) [exp((-1/20).*[0:(T-(21+180))])]]';
signal3.time = [zeros(1,400) [exp((-1/20).*[0:(T-(21+380))])]]';
signal4.time = [zeros(1,600) [exp((-1/20).*[0:(T-(21+580))])]]';
[I,J] = ndgrid(1:m,1:m); r = 5; 
signal1.space = double((I-m/4).^2+(J-m/4).^2<r^2); % what is larger than r^2 is 1 (white), what is smaller is 0 (black)
signal2.space = double((I-3*m/4).^2+(J-3*m/4).^2<r^2); % what is larger than r^2 is 1 (white), what is smaller is 0 (black)
% signal2.space = double((I-m/5).^2+(J-m/5).^2<r^2); % what is larger than r^2 is 1 (white), what is smaller is 0 (black)
signal3.space = double((I-3*m/4).^2+(J-m/4).^2<r^2); % what is larger than r^2 is 1 (white), what is smaller is 0 (black)
signal4.space = double((I-m/4).^2+(J-3*m/4).^2<r^2); % what is larger than r^2 is 1 (white), what is smaller is 0 (black)

figure; suptitle('signal 1 - time and spatial features');
subplot 121
imagesc(signal1.space); colormap(gray);
subplot 122
plot(signal1.time);
figure; suptitle('signal 2 - time and spatial features');
subplot 121
imagesc(signal2.space); colormap(gray);
subplot 122
plot(signal2.time);
figure; suptitle('signal 3 - time and spatial features');
subplot 121
imagesc(signal3.space); colormap(gray);
subplot 122
plot(signal3.time);
figure; suptitle('signal 4 - time and spatial features');
subplot 121
imagesc(signal4.space); colormap(gray);
subplot 122
plot(signal4.time);
%% noise temporal and spatial components
noise1.time = normrnd(0,0.1,T,1);
[I1,I2] = ndgrid([repmat(linspace(0,2*pi,m/2),1,2)]',[repmat(linspace(0,2*pi,m/2),1,2)]);
noise1.space = cos(I1);
figure; suptitle('noise - time and spatial features');
subplot 121
imagesc(noise1.space); colormap(gray);
subplot 122
plot(noise1.time);

noise2.time = normrnd(0,0.1,T,1);
noise2.space = cos(I2);
figure; suptitle('noise - time and spatial features');
subplot 121
imagesc(noise2.space); colormap(gray);
subplot 122
plot(noise2.time);
t = linspace(0,(T-1)./fs,T);
noise3.time = 0.5*cos(2*pi*3.*t')+normrnd(0,0.1,T,1);
[I1,I2] = ndgrid([repmat(linspace(0,2*pi,m/2),1,2)]',1:m);
noise3.space = cos(I1+I2);
figure; suptitle('noise - time and spatial features');
subplot 121
imagesc(noise3.space); colormap(gray);
subplot 122
plot(noise3.time);

%% construct Z and show 9 frames
for i = 1:T
    Z(:,i) = reshape(signal1.time(i)*signal1.space+signal2.time(i)*signal2.space+...
        signal3.time(i)*signal3.space+signal4.time(i)*signal4.space+...
        noise1.time(i)*noise1.space+noise2.time(i)*noise2.space+noise3.time(i)*noise3.space,[],1);
end
figure; suptitle('the recorded signal (Z) at random frames');
ind2show = randi(1000,9,1);
for i=1:9
    subplot(3,3,i)
    imagesc(reshape(Z(:,ind2show(i)),m,[])); colormap(gray);
    title(['Frame #: ' num2str(ind2show(i))]);
end
%% call tsca function
close all;
NoiseGamma = [repmat(-1:0.2:1,3,1)]';
Alpha = zeros(size(NoiseGamma,1),4);
for i=1:size(NoiseGamma,1)
    [ projected,components,D,Alpha(i,:) ] = tscaFunc(Z,signal1,[noise1,noise2,noise3],[1,NoiseGamma(i,:)],1,1);
end

[ projected,components,D,Alpha ] = tscaFunc(Z,signal1,[noise1,noise2,noise3],[1,-0.25,-0.25,-0.25],1,0);
if sum((projected'-signal1.time).^2)>sum(((-1).*projected'-signal1.time).^2) % if inverted
    figure;imagesc(-1.*reshape(components(:,1),m,[])); title('signal spatial component'); colormap(gray);
else
    figure;imagesc(reshape(components(:,1),m,[])); title('signal spatial component'); colormap(gray);
end
figure;imagesc(reshape(components(:,2),m,[])); title('noise spatial component'); colormap(gray);
figure;imagesc(reshape(components(:,3),m,[])); title('noise spatial component'); colormap(gray);
figure;imagesc(reshape(components(:,4),m,[])); title('noise spatial component'); colormap(gray);
figure;plot(abs(diag(D)),'-*'); title('eigenvalues');
%%
noiseNew.time = eye(T)/T; noiseNew.space = [];
noise2New.time = createToeplitz(3,0.1,1,1,T); noise2New.space = [];
[ projected,components,D,Alpha,output ] = tscaFunc(Z,signal1,[noiseNew noise2New],[1 -0.25*ones(1,2)],100,1);
tscaAnalyze(output,3,[],0,T);
