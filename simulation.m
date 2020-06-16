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
noiseSig = [0 0.1 0.8 2 10];
t = linspace(0,(T-1)/fs,T);
for kk = 1:length(noiseSig) % iterate different noise sig
    runSummary = cell(6,50);
    retinotopicMapTSCA = zeros(m*m,3); % preallocate 40x40x3(hsv)
    retinotopicMapTmax = zeros(m*m,3); % preallocate 40x40x3(hsv)
    retinotopicMapAOF = zeros(m*m,3); % preallocate 40x40x3(hsv)
    retinotopicMapCorr = zeros(m*m,3); % preallocate 40x40x3(hsv)
    retinotopicMapGLM = zeros(m*m,3); % preallocate 40x40x3(hsv)
    retinotopicMapNADAV = zeros(m*m,3); % preallocate 40x40x3(hsv)
    for k=1:size(runSummary,2)
        %% construct signals
        [I,J] = ndgrid(1:m,1:m);
        % locs = [0.25 0.25; 0.25 0.75; 0.5 0.5; 0.75 0.25; 0.75 0.75];
        locs = [0.4 0.4; 0.4 0.6; 0.5 0.5; 0.6 0.4; 0.6 0.6]; ind2plot = [1,5,9,13,17];
        DC = 0.02; r = 4;
%         figure;ind2plot = [1,5,9,13,17];
        amp = (1.5-0.5).*rand(1,5)+0.5; % amplitude distributed U~[0.5 1.5]
        for i=1:5
            signals(i).time = amp(i).*(2.5.*normpdf(0:0.1:(T-1)/10,stim_time(i),1)); % normpdf with peak amplitude 1, times random amplitude
            signals(i).space = MinMaxNorm(imgaussfilt(double((I-m*locs(i,1)).^2+(J-m*locs(i,2)).^2<r^2),sqrt(r)))+DC; % what is larger than r^2 is 1 (white), what is smaller is 0 (black)
            %     signals(i).space = double((I-m*locs(i,1)).^2+(J-m*locs(i,2)).^2<r^2); % what is larger than r^2 is 1 (white), what is smaller is 0 (black)
%             subplot(3,6,ind2plot(i))
%             imagesc(signals(i).space); colormap(gray);
%             title(['sig. ' num2str(i) ' - spatial']);
%             subplot(3,6,ind2plot(i)+1)
%             plot(t,signals(i).time,'color','k');
%             title(['sig. ' num2str(i) ' - temporal']); xlabel('time [sec]'); ylabel('Amp. [$\mu$V]');
%             xticks(1:2:9);
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
%         figure; ind2plot = [1,3,5];
%         for i=1:length(noises)
%             subplot(3,2,ind2plot(i))
%             imagesc(noises(i).space); colormap(gray);
%             title(['noise ' num2str(i) ' - spatial']);
%             subplot(3,2,ind2plot(i)+1)
%             plot(t,noises(i).time,'k');
%             title(['noise ' num2str(i) ' - temporal']); xlabel('time [sec]'); ylabel('Amp. [$\mu$V]');
%         end
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
        % for i=1:size(Z,1) % iterate pixels and add exponential trend)
        %     ZZ(i,:) = Z(i,:)+(0.5.*exp(-t./4));
        % end
%         figure; %suptitle('the recorded signal (Z) at random frames');
%         ind2show = sort([100;501;899;randi(1000,7,1)]);
%         for i=1:length(ind2show)
%             subplot(2,5,i)
%             imshow(reshape(Z(:,ind2show(i)),m,[])); colormap(gray);
%             title(['Frame ' num2str(ind2show(i))]);
%         end
        theoreticalSig = zeros(length(signals(1).time),length(signals));
        for i=1:length(signals)
            theoreticalSig(:,i) = normpdf(0:0.1:(T-1)/10,stim_time(i),1);
        end
        [ZZ,ZZZ,ZZZZ,betas] = GLM_VSDI(Z,[0.67 3],theoreticalSig);
%         ind2show = sort([100;501;899;randi(1000,7,1)]);
%         for i=1:length(ind2show)
%             subplot(2,5,i)
%             imshow(reshape(ZZZZ(:,ind2show(i)),m,[])); colormap(gray);
%             title(['Frame ' num2str(ind2show(i))]);
%         end
        
        %% ============================= TSCA ===================================
        noiseNew.time = eye(T)/T; noiseNew.space = []; mapTSCA = zeros(m,m,length(signals));
        noise2New.time = createToeplitz(3,0.1,1,1,T); noise2New.space = [];
        noise3New.time = createToeplitz(0.67,0.1,1,1,T); noise3New.space = [];
        for i=1:length(signals)
            sig.time = theoreticalSig(:,i)';
            [projected,components,D,Alpha,output] = tscaFunc(Z - mean(Z(:)),sig,[noiseNew noise2New noise3New],[1 -0.2*ones(1,3)],100,1);
            %     tscaAnalyze(output,3,[],0,T);
            [~,I] = max(corr(abs(projected(1:4,:)'),theoreticalSig(:,i))); % get the index of component with highest correlation to original signal
            mapTSCA(:,:,i) = MinMaxNorm(abs(reshape(components(:,I),m,m)));
            [tscaMSE(i),tscaPSNR(i),tscaCNR(i),tscaMSSIM(i),tscaCorr(i),tscaCP(i)] = getPerformance(mapTSCA(:,:,i),signals(i).space,signals(i).space>10*DC,signals(i).space<=10*DC);
        end
        retinotopicMapTSCA = retinotopicMapTSCA + retinotopicMapFromIndividualMaps(mapTSCA,0,'TSCA',85);
        clusterEvalTSCA(k) = ClusterEvaluationSim(mapTSCA);
        % [maxmap,maxind] = max(map,[],3); maxmap = maxmap(:); maxind = maxind(:);
        % map22 = colors(maxind,:); map22(maxmap(:)< prctile(maxmap(:),70),:) = repmat([0 0 0],length(find(maxmap(:)< prctile(maxmap(:),70))),1);
        % retinotopicMap = hsv2rgb(maxind/length(unique(maxind)),ones(size(maxind)),(maxmap-min(maxmap))./(max(maxmap)-min(maxmap)));
        % figure;imagesc(reshape(retinotopicMap,40,40,3));
        
        %% calc performance measures between theoretical signal and Z
        for i=1:length(signals)
            [~,I] = max(signals(i).time);
            mapAOF(:,:,i) = MinMaxNorm(reshape(mean(Z(:,I-25:I+25),2),40,40));
            [origMSE(i),origPSNR(i),origCNR(i),origMSSIM(i),origCorr(i),origCP(i)] = getPerformance( mapAOF(:,:,i),signals(i).space,signals(i).space>10*DC,signals(i).space<=10*DC);
            %     figure;imagesc(reshape(mean(Z(:,I-25:I+25),2),40,40));
        end
        retinotopicMapAOF = retinotopicMapAOF + retinotopicMapFromIndividualMaps(mapAOF,0,'AOF',85);
        clusterEvalAOF(k) = ClusterEvaluationSim(mapAOF);
        
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
            if max(temp(:))~=min(temp(:)) % perform maxminnorm if possible
                mapTmax(:,:,i) = MinMaxNorm(reshape(temp,40,40));
            else
                mapTmax(:,:,i) = reshape(temp,40,40);
            end
            [T_MSE(i),T_PSNR(i),T_CNR(i),T_MSSIM(i),T_corr(i),T_CP(i)] = getPerformance(mapTmax(:,:,i),signals(i).space,signals(i).space>10*DC,signals(i).space<=10*DC);
        end
        retinotopicMapTmax = retinotopicMapTmax + retinotopicMapFromIndividualMaps(mapTmax,0,'Tmax',85);
        clusterEvalTmax(k) = ClusterEvaluationSim(mapTmax);
        % tmax2 = Anew';
        % tmax2(tmax2>=0) = (tmax2(tmax2>=0)./100+2)./2; tmax22 = zeros(length(tmax2),3);% transform 0,200,400,600,800 to 1,2,3,4,5
        % tmax22(tmax2>0,:) = colors(tmax2(tmax2>0),:);
        % figure;imagesc(reshape(tmax22,40,40,3));
        
        %% ======================== Correlation method =========================
        normalizedZ=Z;
        mapCorr = normalizedZ*theoreticalSig; mapCorr = rshp(mapCorr);
        for i=1:length(signals)
            mapCorr(:,:,i) = MinMaxNorm(mapCorr(:,:,i));
            [Corr_MSE(i),Corr_PSNR(i),Corr_CNR(i),Corr_MSSIM(i),Corr_corr(i),Corr_CP(i)] = getPerformance(mapCorr(:,:,i),signals(i).space,signals(i).space>10*DC,signals(i).space<=10*DC);
        end
        retinotopicMapCorr = retinotopicMapCorr + retinotopicMapFromIndividualMaps(mapCorr,0,'corr',85);
        clusterEvalCorr(k) = ClusterEvaluationSim(mapCorr);
        
        %%  =========================== GLM method =============================
        mapGLM = rshp(betas(6:end,:)');
        for i=1:length(signals)
            mapGLM(:,:,i) = MinMaxNorm(mapGLM(:,:,i));
            [GLM_MSE(i),GLM_PSNR(i),GLM_CNR(i),GLM_MSSIM(i),GLM_corr(i),GLM_CP(i)] = getPerformance(mapGLM(:,:,i),signals(i).space,signals(i).space>10*DC,signals(i).space<=10*DC);
        end
        retinotopicMapGLM = retinotopicMapGLM + retinotopicMapFromIndividualMaps(mapGLM,0,'glm',85);
        clusterEvalGLM(k) = ClusterEvaluationSim(mapGLM);

        %% ========================== Nadav's method ===========================
        for i=1:length(signals)
            thisSig = Z(:,stim_time(i)*10-99:stim_time(i)*10+100);
            [ind_mx1 w_mx i_rx w_rx]=xplmts(thisSig,[],mean(thisSig(:))+0.2*std(thisSig(:)),[0.9 1.1],[mean(thisSig(:))+1*std(thisSig(:)) 30], 10);
            if max(ind_mx1(:))~=min(ind_mx1(:)) % if can perform min max norm
                mapMPT(:,:,i) = MinMaxNorm(reshape(ind_mx1,40,40));
            else
                mapMPT(:,:,i) = reshape(ind_mx1,40,40);
            end
            %     figure;imagesc(reshape(mapN(:,:,i),40,40));colormap(gray);
            [nadavMSE(i),nadavPSNR(i),nadavCNR(i),nadavMSSIM(i),nadavCorr(i),nadavCP(i)] = getPerformance(mapMPT(:,:,i),signals(i).space,signals(i).space>10*DC,signals(i).space<=10*DC);
        end
        retinotopicMapNADAV = retinotopicMapNADAV + retinotopicMapFromIndividualMaps(mapMPT,0,'nadav',85);
        clusterEvalNadav(k) = ClusterEvaluationSim(mapMPT);
        %% finish up this iteration
        % close all;
        runSummary{1,k} = [origMSE;origPSNR;origCNR;origMSSIM;origCorr;origCP];
        runSummary{2,k} = [nadavMSE;nadavPSNR;nadavCNR;nadavMSSIM;nadavCorr;nadavCP];
        runSummary{3,k} = [tscaMSE;tscaPSNR;tscaCNR;tscaMSSIM;tscaCorr;tscaCP];
        runSummary{4,k} = [T_MSE;T_PSNR;T_CNR;T_MSSIM;T_corr;T_CP];
        runSummary{5,k} = [Corr_MSE;Corr_PSNR;Corr_CNR;Corr_MSSIM;Corr_corr;Corr_CP];
        runSummary{6,k} = [GLM_MSE;GLM_PSNR;GLM_CNR;GLM_MSSIM;GLM_corr;GLM_CP];       
    end
    %% Compare results (all runs)
    ORIG = cat(3,runSummary{1,:});
    NADAV = cat(3,runSummary{2,:}); NADAV(isnan(NADAV)) = 0;
    TSCA = cat(3,runSummary{3,:});
    Tmax = cat(3,runSummary{4,:}); Tmax(isinf(Tmax)) = 100; Tmax(isnan(Tmax)) = 0;
    Corr = cat(3,runSummary{5,:});
    GLM = cat(3,runSummary{6,:});
    clusterEvalAll = {clusterEvalTSCA,clusterEvalAOF,clusterEvalTmax,clusterEvalCorr,clusterEvalGLM,clusterEvalNadav};
    Title = {'MSE','PSNR','CNR','MSSIM','Corr','CP'};
    Title2 = {'TSCA','Tmax','AOF','Corr','GLM','MPT'};
    retMaps = {retinotopicMapTSCA./k,retinotopicMapTmax./k,retinotopicMapAOF./k,retinotopicMapCorr./k,retinotopicMapGLM./k,retinotopicMapNADAV./k};
    figure(kk*10);
    figure(kk*100);
    for i=1:6 % iterate the 6 performance measures
        figure(kk*10);
        subplot(2,3,i)
        boxplot([squeeze(mean(TSCA(i,:,:),2)) squeeze(mean(Tmax(i,:,:),2)) squeeze(mean(ORIG(i,:,:),2))...
            squeeze(mean(Corr(i,:,:),2)) squeeze(mean(GLM(i,:,:),2)) squeeze(mean(NADAV(i,:,:),2))]...
            ,char({'TSCA';'T';'AOF';'Corr';'GLM';'MPT'}));
        title(Title{i});
        figure(kk*100);
        subplot(2,3,i)
        imagesc(reshape(retMaps{i},40,40,3)); title([Title2{i} ' Retinotopic Map']);
    end
    thisSNR_Summary = {TSCA,Tmax,ORIG,Corr,GLM,NADAV,retMaps,clusterEvalAll};
    save(fullfile('C:\Users\Ori\Desktop\Ori\2nd degree\matlab codez\vsdi - matlab\simulation results',...
        ['SimulationSummary_NoiseSig=' num2str(noiseSig(kk)) '.mat']),'thisSNR_Summary');
end