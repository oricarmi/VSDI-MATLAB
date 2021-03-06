%% calculate and show adjusted Si and DBI all maps
path = "C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\comparision results 2";
files = dir(path);
AllRMats = cell(1,7);
AllDBI = cell(1,7);
AllS = cell(1,7);
AllmedS = cell(1,7);
global params brn
for i=1:length(files)-1 % iterate files (except last one which is NIR)  
    if ~isfolder(files(i).name) && ( ( contains(files(i).name,'loc 8') && ~contains(files(i).name,'n=3') ...
            && ~contains(files(i).name,'191218') ) ) || ( contains(files(i).name,'200121') )
    load(fullfile(files(i).folder,files(i).name)); % load summary
    result = Summary2.result;
    fn = fieldnames(result);
    params = Summary2.params;
    RR = mean(params.experiment.optimalMaps.orig,3)>prctile(reshape(mean(Summary2.params.experiment.optimalMaps.orig,3),[],1),85);
    for j=1:length(fn)-1 % iterate the methods and the optimal maps and perform cluster analysis between the maps
        [~,~,maxind{j}] = retinotopicMapFromIndividualMaps(result.(fn{j}).maps,0,'',93);
        maxind{j} = maxind{j}.*RR;
        for k=1:Summary2.params.experiment.N % iterate the maps of this method
            [row,col] = find(maxind{j}==k);
            X{k,j} = [row,col];
        end
    end
    for j=1:length(fn)-1
        [s{j},h] = silhouette(cat(1,X{:,j}),makeClustVector(cellfun(@(x) size(x,1),X(:,j)))');
        NN = histcounts(makeClustVector(cellfun(@(x) size(x,1),X(:,j)))'); % size of each cluster
        expectedClustSize = length(makeClustVector(cellfun(@(x) size(x,1),X(:,j)))')/params.experiment.N; % expected cluster size
        NCM = abs(params.experiment.N-length(unique(NN)))/params.experiment.N; % number of clusters mismatch
        SP = abs(NN - expectedClustSize)/expectedClustSize; % size penalty of each cluster
        TP = SP + NCM; % total penalty
        disp([fn{j} ' ' num2str(mean(TP))])
        TP_allPoints = repelem(TP,NN);
        s{j} = s{j} - TP_allPoints';
%         s{j} = 2*(s{j}-min(s{j}))/(max(s{j})-min(s{j}))-1;
        AllS{j} = [AllS{j};s{j}];
        AllmedS{j} = [AllmedS{j};median(s{j})];
        [R{j},~,~,~,~,DBI{j}] = ClusterSimilarity(X(:,j));
        if contains(files(i).name,'loc 8')
            AllRMats{j} = cat(3,AllRMats{j},R{j});
        end
        DBI{j} = DBI{j} + mean(TP);
        AllDBI{j} = [AllDBI{j};DBI{j}];
    end
%     close all
    clear X maxind R DBI s
%     figure;
%     for i=1:7
%         subplot(3,3,i)
%         imagesc(maxind{i});
%     end
    end
end
figure;
% AllS = cellfun(@(x) max(x, -ones(size(x))),AllS,'uniformoutput',false);
boxplot([AllS{2};AllS{3};AllS{4};AllS{5};AllS{6};AllS{7};AllS{1}],[makeClustVector(cellfun(@(x) size(x,1),AllS))]','labels',char({'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT';'TSCA+GLM'}));
title('Silhouette index');
ylabel('Si');
figure;
boxplot([AllmedS{2};AllmedS{3};AllmedS{4};AllmedS{5};AllmedS{6};AllmedS{7};AllmedS{1}],[makeClustVector(cellfun(@(x) size(x,1),AllmedS))]','labels',char({'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT';'TSCA+GLM'}));
title('Median Silhouette index');
ylabel('median Si');
figure;
boxplot([AllDBI{2},AllDBI{3},AllDBI{4},AllDBI{5},AllDBI{6},AllDBI{7},AllDBI{1}],char({'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT';'TSCA+GLM'}));
title('Davies-Bouldin Index');
ylabel('DB Index');
%%
%% calculate and show adjusted Si and DBI and the retmaps all of a particular map
path = 'E:\comparision results 2';
files = dir(path);
AllRMats = cell(1,7);
AllDBI = cell(1,7);
AllS = cell(1,7);
AllmedS = cell(1,7);
global params brn
retmapsAll = figure("name","retmapsAll");
for i=1:length(files)-1 % iterate files (except last one which is NIR)  
    if ~isfolder(files(i).name) && ( contains(files(i).name,'191119') )
    load(fullfile(files(i).folder,files(i).name)); % load summary
    result = Summary2.result;
    fn = fieldnames(result);
    params = Summary2.params;
    

    
    RR = mean(params.experiment.optimalMaps.orig,3)>prctile(reshape(mean(Summary2.params.experiment.optimalMaps.orig,3),[],1),85);
    for j=1:length(fn)-1 % iterate the methods and the optimal maps and perform cluster analysis between the maps
        [~,~,maxind{j}] = retinotopicMapFromIndividualMaps(result.(fn{j}).maps,0,'',93);
        
        maxind{j} = maxind{j}.*RR;
        for k=1:Summary2.params.experiment.N % iterate the maps of this method
            [row,col] = find(maxind{j}==k);
            X{k,j} = [row,col];
        end
    end
    for j=1:length(fn)-1
        [s{j},h] = silhouette(cat(1,X{:,j}),makeClustVector(cellfun(@(x) size(x,1),X(:,j)))');
        NN = histcounts(makeClustVector(cellfun(@(x) size(x,1),X(:,j)))'); % size of each cluster
        expectedClustSize = length(makeClustVector(cellfun(@(x) size(x,1),X(:,j)))')/params.experiment.N; % expected cluster size
        NCM = (params.experiment.N-length(unique(NN))); % number of clusters mismatch
        SP = abs(NN - expectedClustSize)/expectedClustSize; % size penalty of each cluster
        TP = SP + NCM; % total penalty
        disp([fn{j} ' ' num2str(mean(TP))])
        TP_allPoints = repelem(TP,NN);
        s{j} = s{j} - TP_allPoints';
        AllS{j} = [AllS{j};s{j}];
        AllmedS{j} = [AllmedS{j};median(s{j})];
        [R{j},~,~,~,~,DBI{j}] = ClusterSimilarity(X(:,j));
        if contains(files(i).name,'loc 8')
            AllRMats{j} = cat(3,AllRMats{j},R{j});
        end
        DBI{j} = DBI{j} + mean(TP);
        AllDBI{j} = [AllDBI{j};DBI{j}];
    end % iterate the methods, calculate adjusted silhouette and DBI
%     close all
    clear X maxind R DBI s
%     figure;
%     for i=1:7
%         subplot(3,3,i)
%         imagesc(maxind{i});
%     end
    end
end
figure;
% AllS = cellfun(@(x) max(x, -ones(size(x))),AllS,'uniformoutput',false);
boxplot([AllS{2};AllS{3};AllS{4};AllS{5};AllS{6};AllS{7};AllS{1}],[makeClustVector(cellfun(@(x) size(x,1),AllS))]','labels',char({'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT';'TSCA+GLM'}));
% title('Silhouette index');
ylabel('Si');

figure;
subplot 121
boxplot([AllmedS{1};AllmedS{2};AllmedS{3};AllmedS{4};AllmedS{5};AllmedS{6};AllmedS{7}],[makeClustVector(cellfun(@(x) size(x,1),AllmedS))]','labels',char({'T&G';'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'}),'symbol', '');
% title('Median Silhouette index');
ylabel('median Adjusted Si');
subplot 122
boxplot([AllDBI{1},AllDBI{2},AllDBI{3},AllDBI{4},AllDBI{5},AllDBI{6},AllDBI{7}],char({'T&G';'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'}),'symbol', '');
ylabel('Adjusted DBI');
% title('Davies-Bouldin Index');

for ii=1:length(fn)-1 % iterate the methods    
    figure(retmapsAll); subplot(3,3,ii)
    [~,r] = retinotopicMapFromIndividualMaps(result.(fn{ii}).maps,0,fn{ii},92);
    imf2(r);
    if contains(fn{ii},'nadav','ignorecase',true)
        title('M.P.T');
    else
        title(fn{ii});
    end
end   % plot all retmaps in one figure
%% statistical analysis
vals = [AllDBI{1};AllDBI{2};AllDBI{3};AllDBI{4};AllDBI{5};AllDBI{6};AllDBI{7}];
ttl = ["TSCA+GLM";"TSCA";"Tmax";"AOF";"Corr";"GLM";"MPT"];
groups = convertStringsToChars(repelem(ttl,8));
[p,t,stats] = anova1(vals,groups,'off');
[cdbi,m,h,nms] = multcompare(stats);
vals = [AllmedS{1};AllmedS{2};AllmedS{3};AllmedS{4};AllmedS{5};AllmedS{6};AllmedS{7}];
ttl = ["TSCA+GLM";"TSCA";"Tmax";"AOF";"Corr";"GLM";"MPT"];
groups = convertStringsToChars(repelem(ttl,8));
[p,t,stats] = anova1(vals,groups,'off');
[cs,m,h,nms] = multcompare(stats);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% show simulation results (without running it all over again)
path = 'C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\simulation results 2';
files = dir(path);
figure(112);
for m=6:6
    load(fullfile(files(m).folder,files(m).name)); 
    TSCA = thisSNR_Summary{1};
    Tmax = thisSNR_Summary{2};
    ORIG = thisSNR_Summary{3};
    Corr = thisSNR_Summary{4};
    GLM = thisSNR_Summary{5};
    NADAV = thisSNR_Summary{6};
    TSCAwGLM = thisSNR_Summary{7};
    Title = {'MSE','PSNR','CNR','MSSIM','Pearson''s Corr','CP'};
    Title2 = {'TSCA & GLM','TSCA','Tmax','AOF','Corr','GLM','MPT'};
    Title3 = {'(1)','(2)','(3)','(4)','(5)','(6)','(7)','(8)'};
    retMaps = thisSNR_Summary{8};
    clusterEvalDBI = thisSNR_Summary{9};
    clusterEvalS = thisSNR_Summary{10};
    figure(m*100);
    figure(m*1000);
    for i=1:6 % iterate the 6 performance measures
        figure(m*100);
        subplot(2,3,i)
        boxplot([squeeze(mean(TSCAwGLM(i,:,:),2)) squeeze(mean(TSCA(i,:,:),2)) squeeze(mean(Tmax(i,:,:),2)) squeeze(mean(ORIG(i,:,:),2))...
            squeeze(mean(Corr(i,:,:),2)) squeeze(mean(GLM(i,:,:),2)) squeeze(mean(NADAV(i,:,:),2))]...
            ,char({'T&G';'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'}),'symbol', '');
        ylabel(Title{i});
        box off;
        
%         figure(m*1000);
%         subplot(2,3,i)
%         imagesc(reshape(retMaps{i},40,40,3)); title(Title3{i});%title([Title2{i} ' Retinotopic Map']);
    end
    figure(m*10000);
    for j=1:7 % iterate method
        for k=1:100 % iterate repitition
            if length(clusterEvalS{j}(k).s)~=240
                clusterEvalS{j}(k).s = [clusterEvalS{j}(k).s;-1*ones(240-length(clusterEvalS{j}(k).s),1)];
            end
        end
    end
    boxplot([median([clusterEvalS{7}.s])' median([clusterEvalS{1}.s])' median([clusterEvalS{2}.s])' median([clusterEvalS{3}.s])'...
        median([clusterEvalS{4}.s])' median([clusterEvalS{5}.s])' median([clusterEvalS{6}.s])']...
        ,char({'T&G';'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'}),'symbol', '');
    ylabel('Silhouette Index');
    figure(m*1000+1);
    boxplot([([clusterEvalDBI{7}.DBI])' ([clusterEvalDBI{1}.DBI])' ([clusterEvalDBI{2}.DBI])' ([clusterEvalDBI{3}.DBI])'...
        ([clusterEvalDBI{4}.DBI])' ([clusterEvalDBI{5}.DBI])' ([clusterEvalDBI{6}.DBI])']...
        ,char({'T&G';'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'}),'symbol', '');
    ylabel('DBI');
    
    if m==6
        for kk=2:7
            figure(112);
            subplot(2,4,kk)
            imagesc(reshape(retMaps{kk-1},40,40,3)); title([Title3{kk}]);
            xlabel('pixels');ylabel('pixels');
        end
        subplot(2,4,1)
        imagesc(reshape(retMaps{7},40,40,3)); title([Title3{1}]);
        xlabel('pixels');ylabel('pixels');
        subplot(2,4,8)
        title([Title3{8}]);
    end
    
end
%%
%% show simulation results (without running it all over again) - supplementary material
path = 'C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\simulation results 2';
files = dir(path);
for m=3:length(files)
    load(fullfile(files(m).folder,files(m).name)); 
    TSCA = thisSNR_Summary{1};
    Tmax = thisSNR_Summary{2};
    ORIG = thisSNR_Summary{3};
    Corr = thisSNR_Summary{4};
    GLM = thisSNR_Summary{5};
    NADAV = thisSNR_Summary{6};
    TSCAwGLM = thisSNR_Summary{7};
    Title = {'MSE','PSNR','CNR','MSSIM','Pearson''s Corr','CP'};
    Title2 = {'TSCA & GLM','TSCA','Tmax','AOF','Corr','GLM','MPT'};
    Title3 = {'(1)','(2)','(3)','(4)','(5)','(6)','(7)','(8)'};
    retMaps = thisSNR_Summary{8};
    clusterEvalDBI = thisSNR_Summary{9};
    clusterEvalS = thisSNR_Summary{10};
    figure(m*100);
    figure(m*1000);
    for i=1:6 % iterate the 6 performance measures
        figure(m*100);
        subplot(2,3,i)
        boxplot([squeeze(mean(TSCAwGLM(i,:,:),2)) squeeze(mean(TSCA(i,:,:),2)) squeeze(mean(Tmax(i,:,:),2)) squeeze(mean(ORIG(i,:,:),2))...
            squeeze(mean(Corr(i,:,:),2)) squeeze(mean(GLM(i,:,:),2)) squeeze(mean(NADAV(i,:,:),2))]...
            ,char({'T&G';'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'}),'symbol', '');
        ylabel(Title{i});
        box off;
    end
    
    % fix silhouette indexes
    for j=1:7 % iterate method
        for k=1:100 % iterate repitition
            if length(clusterEvalS{j}(k).s)~=240
                clusterEvalS{j}(k).s = [clusterEvalS{j}(k).s;-1*ones(240-length(clusterEvalS{j}(k).s),1)];
            end
        end
    end
    % end fix

    % iterate methods and plot retinotopic maps 
    for kk=2:7
        figure(m*1000);
        subplot(2,4,kk)
        imagesc(reshape(retMaps{kk-1},40,40,3)); title([Title3{kk}]);
        xlabel('pixels');ylabel('pixels');
    end
    subplot(2,4,1) % plot tsca+glm in first subplot
    imagesc(reshape(retMaps{7},40,40,3)); title([Title3{1}]);
    xlabel('pixels');ylabel('pixels');
    % end
    % boxblot of silhouette index (last subplot, #8)
    set(0,'DefaultTextFontSize',6);
    subplot(2,4,8) 
    boxplot([median([clusterEvalS{7}.s])' median([clusterEvalS{1}.s])' median([clusterEvalS{2}.s])' median([clusterEvalS{3}.s])'...
    median([clusterEvalS{4}.s])' median([clusterEvalS{5}.s])' median([clusterEvalS{6}.s])']...
    ,char({'T&G';'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'}),'symbol', '');
    ylabel('Silhouette Index');title(Title3{8});box off;ylim([-1,1]);
    set(0,'DefaultTextFontSize',14);
    % end
end
%% calculate and show adjusted Si and DBI for supplementary material (optimal maps)
path = 'E:\comparision results 2';
files = dir(path);
AllRMats = cell(1,8);
AllDBI = cell(1,8);
AllS = cell(1,8);
AllmedS = cell(1,8);
global params brn
for i=1:length(files)-1 % iterate files (except last one which is NIR)  
    if ~isfolder(files(i).name) && ( contains(files(i).name,'191119') )
    load(fullfile(files(i).folder,files(i).name)); % load summary
    result = Summary2.result;
    fn = fieldnames(result);
    params = Summary2.params;
    
    
    RR = mean(params.experiment.optimalMaps.orig,3)>prctile(reshape(mean(Summary2.params.experiment.optimalMaps.orig,3),[],1),85);
    for j=1:length(fn) % iterate the methods and the optimal maps and perform cluster analysis between the maps
        if j<length(fn)
            [~,~,maxind{j}] = retinotopicMapFromIndividualMaps(result.(fn{j}).maps,0,'',93);
        else
            [~,~,maxind{j}] = retinotopicMapFromIndividualMaps(params.experiment.optimalMaps.orig,0,'',93);
        end
        maxind{j} = maxind{j}.*RR;
        for k=1:Summary2.params.experiment.N % iterate the maps of this method
            [row,col] = find(maxind{j}==k);
            X{k,j} = [row,col];
        end
    end
    for j=1:length(fn)
        [s{j},h] = silhouette(cat(1,X{:,j}),makeClustVector(cellfun(@(x) size(x,1),X(:,j)))');
        NN = histcounts(makeClustVector(cellfun(@(x) size(x,1),X(:,j)))'); % size of each cluster
        expectedClustSize = length(makeClustVector(cellfun(@(x) size(x,1),X(:,j)))')/params.experiment.N; % expected cluster size
        NCM = (params.experiment.N-length(unique(NN))); % number of clusters mismatch
        SP = abs(NN - expectedClustSize)/expectedClustSize; % size penalty of each cluster
        TP = SP + NCM; % total penalty
        disp([fn{j} ' ' num2str(mean(TP))])
        TP_allPoints = repelem(TP,NN);
        s{j} = s{j} - TP_allPoints';
        AllS{j} = [AllS{j};s{j}];
        AllmedS{j} = [AllmedS{j};median(s{j})];
        [R{j},~,~,~,~,DBI{j}] = ClusterSimilarity(X(:,j));
        if contains(files(i).name,'loc 8')
            AllRMats{j} = cat(3,AllRMats{j},R{j});
        end
        DBI{j} = DBI{j} + mean(TP);
        AllDBI{j} = [AllDBI{j};DBI{j}];
    end % iterate the methods, calculate adjusted silhouette and DBI
%     close all
    clear X maxind R DBI s
%     figure;
%     for i=1:7
%         subplot(3,3,i)
%         imagesc(maxind{i});
%     end
    end
end

figure;
subplot 121
boxplot([AllmedS{1};AllmedS{2};AllmedS{3};AllmedS{4};AllmedS{5};AllmedS{6};AllmedS{7};AllmedS{8}],[makeClustVector(cellfun(@(x) size(x,1),AllmedS))]','labels',char({'T&G';'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT';'Manual'}),'symbol', '');
% title('Median Silhouette index');
ylabel('median Adjusted Si');
box off
subplot 122
boxplot([AllDBI{1},AllDBI{2},AllDBI{3},AllDBI{4},AllDBI{5},AllDBI{6},AllDBI{7},AllDBI{8}],char({'T&G';'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT';'Manual'}),'symbol', '');
ylabel('Adjusted DBI');
box off
% title('Davies-Bouldin Index');
%%
% hi Lamberto once again