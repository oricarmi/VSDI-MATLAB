addpath('C:\Users\orica\Dropbox\fcns_and_decript');
addpath('C:\Users\orica\Dropbox\master degree\codes');
path = 'C:\Users\Ori\Desktop\Ori\2nd degree\matlab codez\vsdi - matlab\comparison results';
path = 'C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\comparison results';
files = dir(path);
global brn lgn ump
brn = zeros(270,327);
for i=7:length(files) % iterate files
    if contains(files(i).name,'181218')
        continue
    end
    load(fullfile(files(i).folder,files(i).name)); % load summary
    result = Summary.result;
    paramz = Summary.params;
    optimalMaps = Summary.params.experiment.optimalMaps.orig;
    N = Summary.params.experiment.N;
    X = cell(N,7); % cell array of clusters (rows - maps, cols - methods)
    fn = fieldnames(result);
    RR = mean(Summary.params.experiment.optimalMaps.orig,3)>prctile(reshape(mean(Summary.params.experiment.optimalMaps.orig,3),[],1),80);
    for j=1:length(fn)+1 % iterate the methods and the optimal maps and perform cluster analysis between the maps
        if j~=length(fn)+1 % if not in optimal maps
            for k=1:N % iterate the maps of this method
                thisMap = result.(fn{j}).maps(:,:,k).*double(RR);
                [row,col] = find(thisMap>prctile(reshape(thisMap,[],1),99));
%                 mapp(:,:,k) = thisMap>prctile(reshape(thisMap,[],1),99);
                X{k,j} = [row,col];
            end
            [R,dbs,dbsI,dbns,dbnsI,DBI] = ClusterSimilarity(X(:,j));
            result.(fn{j}).clusterEval.R = R;
            result.(fn{j}).clusterEval.dbs = dbs;
            result.(fn{j}).clusterEval.dbsI = dbsI;
            result.(fn{j}).clusterEval.dbns = dbns;
            result.(fn{j}).clusterEval.dbnsI = dbnsI;
            result.(fn{j}).clusterEval.DBI = DBI;
        else % it is optimal maps
            for k=1:N % iterate the optimal maps
                [row,col] = find(optimalMaps(:,:,k)>prctile(reshape(optimalMaps(:,:,k),[],1),99));
%                 mapp(:,:,k) = optimalMaps(:,:,k)>prctile(reshape(optimalMaps(:,:,k),[],1),99);
                X{k,j} = [row,col];
            end
            [R,dbs,dbsI,dbns,dbnsI,DBI] = ClusterSimilarity(X(:,j));
            paramz.experiment.optimalMaps.clusterEval.R = R;
            paramz.experiment.optimalMaps.clusterEval.dbs = dbs;
            paramz.experiment.optimalMaps.clusterEval.dbsI = dbsI;
            paramz.experiment.optimalMaps.clusterEval.dbns = dbns;
            paramz.experiment.optimalMaps.clusterEval.dbnsI = dbnsI;
            paramz.experiment.optimalMaps.clusterEval.DBI = DBI;
        end
    end
    for k=1:N % iterate the maps and perform cluster analysis between the methods
        [R,dbs,dbsI,dbns,dbnsI,DBI] = ClusterSimilarity(X(k,:));
        result.clusterEvalAll{k}.R = R;
        result.clusterEvalAll{k}.dbs = dbs;
        result.clusterEvalAll{k}.dbsI = dbsI;
        result.clusterEvalAll{k}.dbns = dbns;
        result.clusterEvalAll{k}.dbnsI = dbnsI;
        result.clusterEvalAll{k}.DBI = DBI;
    end 
    Summary.result = result;
    Summary.params = paramz;
    uisave('Summary',Summary.description)
end

%% calc actual snr (simulation)
ZSig = zeros(m*m,T); % preallocate memory
for i = 1:T
    ZSig(:,i) = reshape( ...
        signals(1).time(i)*signals(1).space+signals(2).time(i)*signals(2).space+...
        signals(3).time(i)*signals(3).space+signals(4).time(i)*signals(4).space+...
        signals(5).time(i)*signals(5).space,[],1);
end
ZNoise = zeros(m*m,T); % preallocate memory
for i = 1:T
    ZNoise(:,i) = reshape( ...
        noises(1).time(i)*noises(1).space+...
        noises(2).time(i)*noises(2).space+...
        noises(3).time(i)*noises(3).space,[],1);
end
ZSig = reshape(ZSig,40,40,[]);
ZNoise = reshape(ZNoise,40,40,[]);
ZSig = ZSig(0.3*m:0.7*m,0.3*m:0.7*m,:);
ZNoise = ZNoise(0.3*m:0.7*m,0.3*m:0.7*m,:);
SNR = 10*log10((ZSig(:)'*ZSig(:))/(ZNoise(:)'*ZNoise(:)));
%% summarize all cluster evaluations real data
addpath('C:\Users\orica\Dropbox\fcns_and_decript');
addpath('C:\Users\orica\Dropbox\master degree\codes');
path = 'C:\Users\Ori\Desktop\Ori\2nd degree\matlab codez\vsdi - matlab\comparison results';
% path = 'C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\comparison results';
files = dir(path);
global brn lgn ump
brn = zeros(270,327);
allPerformance = []; allDBI = [];
%% summarize all performance measures real data
path = 'C:\Users\Ori\Desktop\Ori\2nd degree\matlab codez\vsdi - matlab\comparison results';
files = dir(path);
totalTSCA = []; totalTmax = []; totalAOF = []; totalCorr = []; totalGLM = []; totalNadav = [];
for i=3:length(files) % iterate files
    load(fullfile(files(i).folder,files(i).name)); % load summary
    result = Summary.result;
    totalTSCA = [totalTSCA;result.TSCA.performance];
    totalTmax = [totalTmax;result.Tmax.performance];
    totalAOF = [totalAOF;result.AOF.performance];
    totalCorr = [totalCorr;result.Corr.performance];
    totalGLM = [totalGLM;result.GLM.performance];
    totalNadav = [totalNadav;result.Nadav.performance];
end
Title = {'MSE','PSNR','CNR','MSSIM','Corr','CP'};
Title2 = {'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'};
figure;
for i=1:size(totalAOF,2) % iterate the 6 performance measures
    subplot(2,3,i);
    boxplot([totalTSCA(:,i),totalTmax(:,i),totalAOF(:,i),totalCorr(:,i),totalGLM(:,i),totalNadav(:,i)],...
        char({'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'}));
    title(Title{i});
end
%% summarize all cluster analysis real data
path = 'C:\Users\Ori\Desktop\Ori\2nd degree\matlab codez\vsdi - matlab\comparison results';
files = dir(path);
totalTSCA = []; totalTmax = []; totalAOF = []; totalCorr = []; totalGLM = []; totalNadav = [];
for i=3:length(files) % iterate files
    load(fullfile(files(i).folder,files(i).name)); % load summary
    result = Summary.result;
    for j=1:Summary.params.experiment.N
        totalTSCA = [totalTSCA;Summary.result.clusterEvalAll{j}.R(1,end)];
        totalTmax = [totalTmax;Summary.result.clusterEvalAll{j}.R(2,end)];
        totalAOF = [totalAOF;Summary.result.clusterEvalAll{j}.R(3,end)];
        totalCorr = [totalCorr;Summary.result.clusterEvalAll{j}.R(4,end)];
        totalGLM = [totalGLM;Summary.result.clusterEvalAll{j}.R(5,end)];
        totalNadav = [totalNadav;Summary.result.clusterEvalAll{j}.R(6,end)];
    end
Title = {'MSE','PSNR','CNR','MSSIM','Corr','CP'};
Title2 = {'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'};
figure;
boxplot([totalTSCA,totalTmax,totalAOF,totalCorr,totalGLM,totalNadav],...
    char({'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'}));
title('Cluster Evaluation');
end

%% save video
figure; 
for l=1:size(ZtoWrite,3)
    imagesc(ZtoWrite(:,:,l));colormap('gray');
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);
%% show simulation results (without running it all over again)
path = 'C:\Users\Ori\Desktop\Ori\2nd degree\matlab codez\vsdi - matlab\simulation results';
files = dir(path);
for m=3:length(files)
    load(fullfile(files(m).folder,files(m).name)); 
    TSCA = thisSNR_Summary{1};
    Tmax = thisSNR_Summary{2};
    ORIG = thisSNR_Summary{3};
    Corr = thisSNR_Summary{4};
    GLM = thisSNR_Summary{5};
    NADAV = thisSNR_Summary{6};
    Title = {'MSE','PSNR','CNR','MSSIM','Pearson''s Corr','CP'};
    Title2 = {'TSCA','Tmax','AOF','Corr','GLM','MPT'};
    retMaps = thisSNR_Summary{7};
    clusterEvalAll = thisSNR_Summary{8};
    figure(m*1000);
    figure(m*10000);
    for i=1:6 % iterate the 6 performance measures
        figure(m*1000);
        subplot(2,3,i)
        boxplot([squeeze(mean(TSCA(i,:,:),2)) squeeze(mean(Tmax(i,:,:),2)) squeeze(mean(ORIG(i,:,:),2))...
            squeeze(mean(Corr(i,:,:),2)) squeeze(mean(GLM(i,:,:),2)) squeeze(mean(NADAV(i,:,:),2))]...
            ,char({'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'}));
        title(Title{i});
        figure(m*10000);
        subplot(2,3,i)
        imagesc(reshape(retMaps{i},40,40,3)); title([Title2{i} ' Retinotopic Map']);
    end
end

%% computation time
addpath("C:\Users\orica\Dropbox\fcns_and_decript");
addpath('C:\Users\orica\Dropbox\master degree\codes');
global fs 
T = 1000; fs = 100;
P = 1600;
m = sqrt(P);
signals = struct;
noises = struct;
stim_time = 10:20:90;
noiseSig = [0 0.1 0.8 2 10];
t = linspace(0,(T-1)/fs,T);
timeTSCA = zeros(1,100);timeTmax = zeros(1,100);timeAOF = zeros(1,100);
timeCorr = zeros(1,100);timeGLM = zeros(1,100);timeMPT = zeros(1,100);
for k=1:100
    %% construct signals
    [I,J] = ndgrid(1:m,1:m);
    locs = [0.4 0.4; 0.4 0.6; 0.5 0.5; 0.6 0.4; 0.6 0.6]; ind2plot = [1,5,9,13,17];
    DC = 0.02; r = 4;
    amp = (1.5-0.5).*rand(1,5)+0.5; % amplitude distributed U~[0.5 1.5]
    for i=1:5
        signals(i).time = amp(i).*(2.5.*normpdf(0:0.1:(T-1)/10,stim_time(i),1)); % normpdf with peak amplitude 1, times random amplitude
        signals(i).space = MinMaxNorm(imgaussfilt(double((I-m*locs(i,1)).^2+(J-m*locs(i,2)).^2<r^2),sqrt(r)))+DC; % what is larger than r^2 is 1 (white), what is smaller is 0 (black)
    end
    %% construct noise
    [I1,I2] = ndgrid([repmat(linspace(0,2*pi,m/2),1,2)]',[repmat(linspace(0,2*pi,m/2),1,2)]);
    noises(1).time = normrnd(0,noiseSig(kk),T,1);
    noises(1).space = MinMaxNorm(normrnd(0,1,m,m))+DC;%cos(I1);
    freqs = 2*pi.*[2.8:0.05:3.2]'; % bandwidth around 3 hz
    phase = 2*pi.*rand(length(freqs),1); % random phase
    amp = noiseSig(kk)*4.*rand(size(phase)); % random amplitude
    noises(2).time = mean(repmat(amp,1,T).*cos(freqs.*t+phase))+normrnd(0,noiseSig(kk)*0.05,1,T); % create almost periodic signal around 3hz
    noises(2).space = MinMaxNorm(cos(3.*I1+2*pi*rand(1)))+DC;
    freqs = 2*pi.*[0.47:0.05:0.87]'; % bandwidth around 0.67 hz
    phase = 2*pi.*rand(length(freqs),1);% random phase
    amp = noiseSig(kk)*4.*rand(size(phase));% random amplitude
    noises(3).time = mean(repmat(amp,1,T).*cos(freqs.*t+phase))+normrnd(0,noiseSig(kk)*0.05,1,T); % create almost periodic signal around 0.67 hz %noiseSig(kk)*cos(2*pi*0.67.*t'+2*pi*rand(1))+normrnd(0,noiseSig(kk)/2,T,1);
    noises(3).space = MinMaxNorm(normrnd(0,1,m,m))+DC;
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
        theoreticalSig = zeros(length(signals(1).time),length(signals));
        for i=1:length(signals)
            theoreticalSig(:,i) = normpdf(0:0.1:(T-1)/10,stim_time(i),1);
        end
%         [ZZ,ZZZ,ZZZZ,betas] = GLM_VSDI(Z,[0.67 3],theoreticalSig);

        
        %% ============================= TSCA ===================================
        tic
        noiseNew.time = eye(T)/T; noiseNew.space = []; mapTSCA = zeros(m,m,length(signals));
        noise2New.time = createToeplitz(3,0.1,1,1,T); noise2New.space = [];
        noise3New.time = createToeplitz(0.67,0.1,1,1,T); noise3New.space = [];
        for i=1:length(signals)
            sig.time = theoreticalSig(:,i)';
            [projected,components,D,Alpha,output] = tscaFunc(Z - mean(Z(:)),sig,[noiseNew noise2New noise3New],[1 -0.2*ones(1,3)],100,1);
            %     tscaAnalyze(output,3,[],0,T);
            [~,I] = max(corr(abs(projected(1:4,:)'),theoreticalSig(:,i))); % get the index of component with highest correlation to original signal
            mapTSCA(:,:,i) = MinMaxNorm(abs(reshape(components(:,I),m,m)));
        end    
        thisTOC = toc;
        timeTSCA(k) = thisTOC;
        %% calc performance measures between theoretical signal and Z
        tic
        for i=1:length(signals)
            [~,I] = max(signals(i).time);
            mapAOF(:,:,i) = MinMaxNorm(reshape(mean(Z(:,I-25:I+25),2),40,40));
        end
        thisTOC = toc;
        timeAOF(k) = thisTOC;
        %% ========================== T_max method =============================
        tic
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
            if max(temp(:))~=min(temp(:)) % perform maxminnorm if possible
                mapTmax(:,:,i) = MinMaxNorm(reshape(temp,40,40));
            else
                mapTmax(:,:,i) = reshape(temp,40,40);
            end
        end 
        thisTOC = toc;
        timeTmax(k) = thisTOC;
        %% ======================== Correlation method =========================
        tic
        mapCorr = Z*theoreticalSig; mapCorr = rshp(mapCorr);
        for i=1:length(signals)
            mapCorr(:,:,i) = MinMaxNorm(mapCorr(:,:,i));
        end
        thisTOC = toc;
        timeCorr(k) = thisTOC;

        %%  =========================== GLM method =============================
        tic
        [ZZ,ZZZ,ZZZZ,betas] = GLM_VSDI(Z,[0.67 3],theoreticalSig);
        mapGLM = rshp(betas(6:end,:)');
        for i=1:length(signals)
            mapGLM(:,:,i) = MinMaxNorm(mapGLM(:,:,i));
        end
        thisTOC = toc;
        timeGLM(k) = thisTOC;
        %% ========================== Nadav's method ===========================
        tic
        for i=1:length(signals)
            thisSig = Z(:,stim_time(i)*10-99:stim_time(i)*10+100);
            [ind_mx1 w_mx i_rx w_rx]=xplmts(thisSig,[],mean(thisSig(:))+0.2*std(thisSig(:)),[0.9 1.1],[mean(thisSig(:))+1*std(thisSig(:)) 30], 10);
            if max(ind_mx1(:))~=min(ind_mx1(:)) % if can perform min max norm
                mapMPT(:,:,i) = MinMaxNorm(reshape(ind_mx1,40,40));
            else
                mapMPT(:,:,i) = reshape(ind_mx1,40,40);
            end
        end
        thisTOC = toc;
        timeMPT(k) = thisTOC;
        %% finish up this iteration
end
%% Run everything again with preprocess normalization being none
addpath("C:\Users\orica\Dropbox\fcns_and_decript");
addpath('C:\Users\orica\Dropbox\master degree\codes');
addpath("C:\Users\Ori\Dropbox\fcns_and_decript")
addpath("C:\Users\Ori\Dropbox\master degree\codes")
% path = 'C:\Users\Ori\Desktop\Ori\2nd degree\matlab codez\vsdi - matlab\comparison results';
path = 'C:\Users\Ori\Desktop\Ori\2nd degree\matlab codez\vsdi - matlab\comparison results';
files = dir(path);
global bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0  basis params
for i=11:length(files)
    close all
    load(fullfile(files(i).folder,files(i).name)); % load summary
    locs = strfind(files(i).name,'_');locs(1) = locs(1)+1; locs(2)=locs(2)-1;
    fname = ['E:\' files(i).name(locs(1):locs(2)) '\m' files(i).name(locs(1):locs(2))];
    n = str2num(files(i).name(locs(2)+2)); 
    if strcmp(files(i).name(locs(1):locs(2)),'200121')
        What = 92;
    elseif strcmp(files(i).name(locs(1):locs(2)),'181218') || strcmp(files(i).name(locs(1):locs(2)),'191119') || strcmp(files(i).name(locs(1):locs(2)),'191218')
        What = 8;
    elseif contains(Summary.description,'Hz','IgnoreCase',true) 
        What = str2num(Summary.description(5))*10+2; 
    else
        What = str2num(Summary.description(5));
    end
    result = struct('TSCAwGLM',struct,'TSCAnoGLM',struct,'Tmax',struct,'AOF',struct,'Corr',struct,'GLM',struct,'Nadav',struct);
    [result.TSCAwGLM.maps,result.Tmax.maps,result.AOF.maps,result.Corr.maps,result.GLM.maps,result.Nadav.maps,tscaBoth] = GenerateMaps(fname,n,What); % 3rd parameter: what, 8=loc8,9=loc9,92=mvngbars2hz
    result.TSCAnoGLM.maps = tscaBoth.noGLM;
    % <----- generate retinotopic maps from the individual maps
    fn = fieldnames(result);
    for j=1:length(fn) % iterate the methods 
        [~,result.(fn{j}).retinotopicMap] = retinotopicMapFromIndividualMaps(result.(fn{j}).maps,0,fn{j});
    end
    % ---->
    params.experiment.optimalMaps = Summary.params.experiment.optimalMaps;
    [result.TSCAwGLM.performance,result.TSCAnoGLM.performance,result.Tmax.performance,result.AOF.performance,result.Corr.performance,result.GLM.performance,result.Nadav.performance] = performanceRealData(result);
    result = clusterEvaluation(result);
    Summary2 = struct('params',params,'result',result,'description',Summary.description); 
    save(['C:\Users\Ori\Desktop\Ori\2nd degree\matlab codez\vsdi - matlab\comparision results 2\' Summary.description],'Summary2')
end
