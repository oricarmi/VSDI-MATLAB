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
    Title2 = {'TSCA and GLM','TSCA','Tmax','AOF','Corr','GLM','MPT'};
    retMaps = thisSNR_Summary{8};
    clusterEvalDBI = thisSNR_Summary{9};
    clusterEvalS = thisSNR_Summary{10};
    figure(m*100);
    for i=1:6 % iterate the 6 performance measures
        figure(m*100);
        subplot(2,3,i)
        boxplot([squeeze(mean(TSCAwGLM(i,:,:),2)) squeeze(mean(TSCA(i,:,:),2)) squeeze(mean(Tmax(i,:,:),2)) squeeze(mean(ORIG(i,:,:),2))...
            squeeze(mean(Corr(i,:,:),2)) squeeze(mean(GLM(i,:,:),2)) squeeze(mean(NADAV(i,:,:),2))]...
            ,char({'T&G';'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'}),'symbol', '');set(gca,'XTickLabelRotation',60);
        ylabel(Title{i});
        box off;
        
%         figure(m*1000);
%         subplot(2,3,i)
%         imagesc(reshape(retMaps{i},40,40,3)); title(Title3{i});%title([Title2{i} ' Retinotopic Map']);
    end
    figure(m*10000);
    subplot 121
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
    ylabel('Silhouette Index');box off; set(gca,'XTickLabelRotation',45);
    subplot 122
    boxplot([([clusterEvalDBI{7}.DBI])' ([clusterEvalDBI{1}.DBI])' ([clusterEvalDBI{2}.DBI])' ([clusterEvalDBI{3}.DBI])'...
        ([clusterEvalDBI{4}.DBI])' ([clusterEvalDBI{5}.DBI])' ([clusterEvalDBI{6}.DBI])']...
        ,char({'T&G';'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'}),'symbol', '');
    ylabel('DBI');box off; set(gca,'XTickLabelRotation',45);
    
    figure(112);
    if m==6
        subplot(2,4,1)
        imagesc(reshape(retMaps{7},40,40,3)); title([Title2{1}]);
        xlabel('pixels');ylabel('pixels');set(gca,'xticklabel',[]);set(gca,'xtick',[]);set(gca,'yticklabel',[]);set(gca,'ytick',[]);
        for kk=2:7
            subplot(2,4,kk)
            imagesc(reshape(retMaps{kk-1},40,40,3)); title([Title2{kk}]);
            xlabel('pixels');ylabel('pixels');
            set(gca,'xticklabel',[]);set(gca,'xtick',[]);set(gca,'yticklabel',[]);set(gca,'ytick',[]);
        end
        subplot(2,4,8);
        imagesc(rshp(retinotopicMap));xlabel('pixels'); ylabel('pixels');title('Expected Map');
        set(gca,'xticklabel',[]);set(gca,'xtick',[]);set(gca,'yticklabel',[]);set(gca,'ytick',[]);
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
    Title2 = {'TSCA and GLM','TSCA','Tmax','AOF','Corr','GLM','MPT'};
    Title3 = {'(1)','(2)','(3)','(4)','(5)','(6)','(7)','(8)'};
    retMaps = thisSNR_Summary{8};
    clusterEvalDBI = thisSNR_Summary{9};
    clusterEvalS = thisSNR_Summary{10};
    figure("name",num2str(m*100),'Position', [100 100 1800 600])
    for i=1:6 % iterate the 6 performance measures
        subplot(2,3,i)
        boxplot([squeeze(mean(TSCAwGLM(i,:,:),2)) squeeze(mean(TSCA(i,:,:),2)) squeeze(mean(Tmax(i,:,:),2)) squeeze(mean(ORIG(i,:,:),2))...
            squeeze(mean(Corr(i,:,:),2)) squeeze(mean(GLM(i,:,:),2)) squeeze(mean(NADAV(i,:,:),2))]...
            ,char({'T&G';'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'}),'symbol', '');
        ylabel(Title{i});
        box off; set(gca,'XTickLabelRotation',45);
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
    figure("name",num2str(m*1000),'Position', [100 100 1800 900])
    for kk=2:7
        subplot(2,4,kk)
        imagesc(reshape(retMaps{kk-1},40,40,3)); title([Title2{kk}]);
        xlabel('pixels');ylabel('pixels');
    end
    subplot(2,4,1) % plot tsca+glm in first subplot
    imagesc(reshape(retMaps{7},40,40,3)); title([Title2{1}]);
    xlabel('pixels');ylabel('pixels');
    subplot(2,4,8);
    imagesc(rshp(retinotopicMap));xlabel('pixels'); ylabel('pixels');title('Expected Map');
    set(gca,'xticklabel',[]);set(gca,'xtick',[]);set(gca,'yticklabel',[]);set(gca,'ytick',[]);
    % end
    % boxblot of silhouette index (last subplot, #8)
%     set(0,'DefaultTextFontSize',6);
    figure("name",num2str(m*10),'Position', [100 100 1600 600])
    subplot 121
    boxplot([median([clusterEvalS{7}.s])' median([clusterEvalS{1}.s])' median([clusterEvalS{2}.s])' median([clusterEvalS{3}.s])'...
    median([clusterEvalS{4}.s])' median([clusterEvalS{5}.s])' median([clusterEvalS{6}.s])']...
    ,char({'T&G';'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'}),'symbol', '');
    ylabel('Silhouette Index');box off;ylim([-1,1]);
    set(gca,'XTickLabelRotation',45);
    subplot 122
    boxplot([([clusterEvalDBI{7}.DBI])' ([clusterEvalDBI{1}.DBI])' ([clusterEvalDBI{2}.DBI])' ([clusterEvalDBI{3}.DBI])'...
        ([clusterEvalDBI{4}.DBI])' ([clusterEvalDBI{5}.DBI])' ([clusterEvalDBI{6}.DBI])']...
        ,char({'T&G';'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'}),'symbol', '');
    ylabel('DBI');box off; set(gca,'XTickLabelRotation',45);
    % end
end
%% show simulation results adjusted dbi and si (without running it all over again) 2 - supplementary material 
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
    % calculate adjusted silhouette and DBI
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


    % iterate methods and plot retinotopic maps 
    for kk=2:7
        figure(m*1000);
        subplot(2,4,kk)
        imagesc(reshape(retMaps{kk-1},40,40,3)); title([Title2{kk}]);
        xlabel('pixels');ylabel('pixels');
    end
    subplot(2,4,1) % plot tsca+glm in first subplot
    imagesc(reshape(retMaps{7},40,40,3)); title([Title2{1}]);
    xlabel('pixels');ylabel('pixels');
    % end
    % boxblot of silhouette index (last subplot, #8)
%     set(0,'DefaultTextFontSize',6);
    figure(m*10);
    subplot 121
    boxplot([median([clusterEvalS{7}.s])' median([clusterEvalS{1}.s])' median([clusterEvalS{2}.s])' median([clusterEvalS{3}.s])'...
    median([clusterEvalS{4}.s])' median([clusterEvalS{5}.s])' median([clusterEvalS{6}.s])']...
    ,char({'T&G';'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'}),'symbol', '');
    ylabel('Silhouette Index');box off;ylim([-1,1]);
    set(gca,'XTickLabelRotation',45);
    subplot 122
    boxplot([([clusterEvalDBI{7}.DBI])' ([clusterEvalDBI{1}.DBI])' ([clusterEvalDBI{2}.DBI])' ([clusterEvalDBI{3}.DBI])'...
        ([clusterEvalDBI{4}.DBI])' ([clusterEvalDBI{5}.DBI])' ([clusterEvalDBI{6}.DBI])']...
        ,char({'T&G';'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'}),'symbol', '');
    ylabel('DBI');box off; set(gca,'XTickLabelRotation',45);
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
%% Plot retmaps of a loaded Summary2 manual maps
result = Summary2.result;
params = Summary2.params;
for i=1:8
    params.experiment.optimalMaps.orig(:,:,i) = MinMaxNorm(params.experiment.optimalMaps.orig(:,:,i));
end
fn = fieldnames(result);
retmapsAll = figure("name","retmapsAll",'Position', [100 0 1400 500]);
% indmapsAll = figure("name","indmapsAll");
% Titles = {'(1)','(2)','(3)','(4)','(5)','(6)','(7)','(8)'};
Titles ={'TSCA and GLM','TSCA','Tmax','AOF','Corr','GLM','MPT','Manual Maps'};
flag=0;
for i=1:length(fn)-1 % iterate the methods  
%     if i==1
%         for iii=1:size(result.(fn{i}).maps,3)
%             result.(fn{i}).maps(:,:,iii) = result.(fn{i}).maps(:,:,iii).*R;
%         end   
%     end
    [~,r] = retinotopicMapFromIndividualMaps(result.(fn{i}).maps,0,fn{i},92);
    figure(retmapsAll); subplot(2,4,i)
    imf2(r);xlabel('');set(gca,'xticklabel',[]);set(gca,'xtick',[]);ylabel('');set(gca,'yticklabel',[]);set(gca,'ytick',[]);
%     imagesc(rshp(result.(fn{i}).retinotopicMap));
    title(Titles{i}); hold on;
    if ~flag
        rectangle('Position',[4.5,5,1,0.1],'FaceColor',[1 1 1],'EdgeColor',[1 1 1]);
        text(4.4,4.8,'1mm','color',[1 1 1],'fontsize',12);
        flag=1;
    end
%     annotation(retmapsAll,'doublearrow',[0.1 0.1],[.1 0]);
%     if contains(fn{i},'nadav','ignorecase',true)
%         title('M.P.T');
%     elseif contains(fn{i},'tscano','ignorecase',true)
%         title('TSCA');
%     elseif contains(fn{i},'tscaw','ignorecase',true)
%         title('TSCA & GLM');
%     else
%         title(fn{i});
%     end
%     r = plotMaps(result.(fn{i}).maps,fn{i},1);
end
subplot(2,4,8);
[~,r] = retinotopicMapFromIndividualMaps(params.experiment.optimalMaps.orig,0,fn{i},98);
imf2(r);xlabel('');set(gca,'xticklabel',[]);set(gca,'xtick',[]);ylabel('');set(gca,'yticklabel',[]);set(gca,'ytick',[]);
%     imagesc(rshp(result.(fn{i}).retinotopicMap));
title(Titles{end}); hold on;
%% supp figure 6b (image restoration metrics between all methods compared to optimal maps)
load('C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\comparision results 2\comparison to manual maps 181218 n=2 for supp material.mat');
figure;
Title = {'MSE','PSNR','CNR','MSSIM','Pearson''s Corr','CP'};
for j=1:6 % iterate six metrics
    subplot(2,3,j);
    boxplot([comparison2manualMaps181218_2{1}(:,j),comparison2manualMaps181218_2{2}(:,j),comparison2manualMaps181218_2{3}(:,j), ...
        comparison2manualMaps181218_2{4}(:,j),comparison2manualMaps181218_2{5}(:,j),comparison2manualMaps181218_2{6}(:,j), ...
        comparison2manualMaps181218_2{7}(:,j)],char({'T&G';'TSCA';'Tmax';'AOF';'Corr';'GLM';'MPT'}),'symbol', '');
    set(gca,'XTickLabelRotation',60);ylabel(Title{j});box off;
end
%% supp figure 2 - dbi and silhouette with decreassing SNR
%% show simulation results (without running it all over again) - supplementary material
path = 'C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\simulation results 2';
files = dir(path);
SNRorder = [5,3,4,7,6];
clusterEvalDBI = zeros(100,5,7);
clusterEvalS = zeros(100,5,7);
for m=1:length(SNRorder)
    load(fullfile(files(SNRorder(m)).folder,files(SNRorder(m)).name)); 
    thisclusterEvalDBI = thisSNR_Summary{9};
    thisclusterEvalS = thisSNR_Summary{10};
    % fix silhouette indexes
    for j=1:7 % iterate method
        for k=1:100 % iterate repitition
            if length(thisclusterEvalS{j}(k).s)~=240
                thisclusterEvalS{j}(k).s = [thisclusterEvalS{j}(k).s;-1*ones(240-length(thisclusterEvalS{j}(k).s),1)];
            end
        end
        clusterEvalDBI(:,m,j) = [thisclusterEvalDBI{j}.DBI]';
        clusterEvalS(:,m,j) = [mean([thisclusterEvalS{j}.s])]';
    end
    % end fix
end
Title2 = {'TSCA and GLM','TSCA','Tmax','AOF','Corr','GLM','MPT'};
for i=1:7
    figure("name",Title2{i},'Position', [100 100 1400 600]);
    suptitle(Title2{i});
    subplot 121;
    errorbar(1:5,nanmean(clusterEvalS2(:,:,i)),nanstd(clusterEvalS2(:,:,i)));
    xlim([0,6]);ylabel('adjusted Silhouette Index');
    xticks(1:5);
    xticklabels({'SNR=inf','SNR=5[dB]','SNR=0[dB]','SNR=-5[dB]','SNR=-10[dB]'});
    set(gca,'XTickLabelRotation',45);
    ylim([-1,1]);
    subplot 122;
    errorbar(1:5,nanmean(clusterEvalDBI2(:,:,i)),nanstd(clusterEvalDBI2(:,:,i)));
    xlim([0,6]);ylabel('adjusted DBI');
    xticks(1:5);
    xticklabels({'SNR=inf','SNR=5[dB]','SNR=0[dB]','SNR=-5[dB]','SNR=-10[dB]'});
    set(gca,'XTickLabelRotation',45);
    ylim([-1,7]);
end
%% Figure 4 main
Title2 = {'TSCA and GLM','TSCA','Tmax','AOF','Corr','GLM','MPT'};
figure("name","All",'Position', [50 100 1450 620]);
for i=1:7
    subplot 121;
    plot(1:5,nanmean(clusterEvalS2(:,:,i)),'-o');hold on;
    subplot 122;
    plot(1:5,nanmean(clusterEvalDBI2(:,:,i)),'-o'); hold on;

end
subplot 121;
xlim([0,6]);ylabel('adjusted Silhouette Index');
xticks(1:5);
xticklabels({'SNR=inf','SNR=5[dB]','SNR=0[dB]','SNR=-5[dB]','SNR=-10[dB]'});
set(gca,'XTickLabelRotation',45);
legend(Title2,'FontSize',12,'location','southwest');box off;
subplot 122;
xlim([0,6]);ylabel('adjusted DBI');
xticks(1:5);
xticklabels({'SNR=inf','SNR=5[dB]','SNR=0[dB]','SNR=-5[dB]','SNR=-10[dB]'});
set(gca,'XTickLabelRotation',45);
legend(Title2,'FontSize',12,'location','northwest');box off;
