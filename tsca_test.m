%% time is average of each frame
noise.time = [mean(cfn{1})]';
noise2.time = [mean(cfn{2})]';
figure;
img = cell(8,1);
for i=3:length(cfn)
    disp(i);
    signal.time = [mean(cfn{i})]';
    Z = cfn{i};
    [ projected,components,D ] = tscaFunc(Z,signal,[noise noise2],[1,0,0],1,1);
    subplot(3,3,i-1);
    img{i-1} = imgaussfilt(reshape(components(:,1),size(brn,1),[]),3.5);
%     if sum((projected'-signal.time).^2)>sum(((-1).*projected'-signal.time).^2) % inverted time approximation, so invert
%         img{i-1} = img{i-1}.*(-1);
%     end
    img{i-1} = sign(projected*signal.time).*img{i-1}; % multipy by sign of inner product to fix inversion
    imagesc(img{i-1}); title(['spatial component - ' lgn(i,:)]); colormap(jet);
%     figure;plot(abs(diag(D)),'-*'); title('eigenvalues');
end
saveas(gcf,'tscatest_BSLN_ARFCT_ALLBRAIN.jpg');
save('tsca_output_new','img');
%% time is maximum of each frame (bad results)
noise.time = [max(cfn{1})]';
figure(1); figure(2);
img = cell(9,1);
for i=2:length(cfn)
    disp(i);
    signal.time = [max(cfn{i})]';
    Z = cfn{i};
    [ projected,components,D ] = tscaFunc(Z,signal,noise,[1,0],1);
    subplot(3,3,i-1);
    img{i-1} = imgaussfilt(reshape(components(:,1),size(brn,1),[]),1.5);
    if sum((projected'-signal.time).^2)>sum(((-1).*projected'-signal.time).^2) % inverted time approximation, so invert
        img{i-1} = img{i-1}.*(-1);
    end
    imagesc(img{i-1}); title(['spatial component - ' lgn(i,:)]); %colormap(gray);
%     figure;plot(abs(diag(D)),'-*'); title('eigenvalues');
end
saveas(gcf,'tscatest_max.jpg');
save('tsca_output_new','img');
%% time is average of each frame (PCA to reduce dimentionality)
noise.time = [mean(cfn{1})]';
figure;
img = cell(8,1);
for i=2:length(cfn)
    disp(i);
    signal.time = [mean(cfn{i})]';
    Z = cfn{i};
    [ projected,components,D ] = tscaFunc(Z,signal,noise,[1,-1],1,1);
    subplot(3,3,i-1);
    img{i-1} = imgaussfilt(reshape(components(:,1),size(c0,1),[]),3);
    if sum((projected'-signal.time).^2)>sum(((-1).*projected'-signal.time).^2) % inverted time approximation, so invert
        img{i-1} = img{i-1}.*(-1);
    end
%     imagesc(img{i-1}); title(['spatial component - ' lgn(i,:)]); colormap(jet);
%     figure;plot(abs(diag(D)),'-*'); title('eigenvalues');
end
saveas(gcf,'tscatest2.jpg');
save('tsca_output_new','img');
%% plot results
mxc_Ori = zeros(numel(img{3}),length(img)-1);
% index2plot = [5,7,3,8,2,9,1,6,4]; % 9 locs
index2plot = [6,2,8,4,5,7,1,5]; % 8 locs
figure;
for i=2:length(img)
%     if i~=4 && i~=5
%         subplot(3,3,index2plot(i-1));
%     else
%         subplot(3,2,index2plot(i-1));
%     end
    switch i
        case 2
            subplot(3,3,i)
  
    end
    imagesc(img{i}); title(['spatial component - ' lgn(i+1,:)]); colormap(jet);
    mxc_Ori(:,i-1) = reshape(img{i},[],1);
end
saveas(gcf,'tscatest_BSLN.jpg');
pc=0; [ry cx mx]=cmbn9(mxc_Ori); %% combine 9 conditions %%
[ds rf cf L]=cmbn84(mxc_Ori, rxc(:,2:end-1), 2); % combine 8 conditions
%% remove max latency not between ranges (time is average of each frame before removal)
noise.time = [mean(cfn{1})]';
% noise2.time = [mean(cfn{2})]';
figure;
img = cell(8,1);
for i=2:length(cfn)
    disp(i);
    signal.time = [mean(cfn{i})]';
%     signal.time = zeros(1,100);
%     signal.time(98:100) = 1;
%     signal.time = repmat(signal.time,1,10);
    Z = cfn{i};
    [ind_mx1]=xplmts(Z, [], 0, [0.1 0.4], [3 10], 0);
%     Z(ind_mx1==0,:) = 0;
    [ projected,components,D ] = tscaFunc(Z,signal,[noise],[1,-1],2,1);
    estimatedSignal = components(:,1);
    if ~isempty(ind_mx1)
        estimatedSignal(ind_mx1==0) = 0;
    end
%     estimatedSignal(ind_mx1==0) = 0;
    img{i-1} = imgaussfilt(reshape(estimatedSignal,size(brn,1),[]),3.5);
%     if sum((projected'-signal.time).^2)>sum(((-1).*projected'-signal.time).^2) % inverted time approximation, so invert
%         img{i-1} = img{i-1}.*(-1);
%     end
    img{i-1} = sign(projected(1,:)*signal.time).*img{i-1}; % multipy by sign of inner product to fix inversion
    imagesc(img{i-1}); title(['spatial component - ' lgn(i,:)]); colormap(jet);
%     figure;plot(abs(diag(D)),'-*'); title('eigenvalues');
end
saveas(gcf,'tscatest_new.jpg');
save('tsca_output_new','img');
%%
c1 = reshape(c0,[],1000);
signal = mean(c1);
% figure;
hold on;
plot(signal);
cfn{1} = normrnd(mean(mean(c1(:,1:100))),std(std(c1(:,1:100))),size(c1));
cfn{2} = c1;
legend('1','2','3');
%%
NoiseGamma = [2.*rand(10,2)-1];
for i=1:size(NoiseGamma,1)
    [status,img,fitted2DG] = tsca_loc8(cfn,1,[0.1 0.3],[2 10],[1 NoiseGamma(i,:)],3,1,3,0);
end
% [status,img] = tsca_loc9(cfn,[0.1 0.3],[3 10],[1 -0.4],1,3,0);
[status,img] = tsca_loc9(cfn,[0.1 0.3],[3 10],[1 -1],3,1,3,0);
[status,img] = tsca_ff(cfn,[0.1 0.2],[2 15],[1 -1],2,1,3.5,0);
%%
%% try to get blood vessels from TSCA
T = length(cfn)*100;
noise(1).time = eye(T)*fs; % autocorrelation matrix of white noise
% noise(2).time = createToeplitz(0.98,0.1,3,[1 1 0.1],T);
% noise(2).time = ones(100)*fs; % autocorrelation matrix of DC
f0 = 6.7; %13.18; %0.68 % 12.89 % 8.3
signal1.time = createToeplitz(f0,0.5,3,[1 0.1 0.1],T);
% Means = cellfun(@(x) mean(x,2),cfn0,'uniformoutput',false);
Z = [];
for i=1:length(cfn)
%     if i==2 
%         continue;
%     end
    Z = [Z cfn0{i}];
end
Z = Z - mean(Z,2);
% Z = cfn0{1}-Means{1}; 
t= linspace(0,length(cfn),100); 
gammas = [1 -0.05];
[projected,components,D,Alphas,tscaStruct] = tscaFunc(Z,signal1,noise,gammas,9,1);
% figure;suptitle(['f0 = ' num2str(f0)]);imagesc(reshape(components(:,end),size(brn0,1),[]));colormap(gray);
tscaAnalyze(tscaStruct,3,['f0 = ' num2str(f0)],1,T);

% figure;suptitle([Title ' - eigen-values']); eigenVals = sort(diag(D),'desc');
% plot(eigenVals(1:100),'*-');xlabel('component #');ylabel('eigen value');
