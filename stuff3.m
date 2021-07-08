%% main figure 9
folder = 'C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\comparision results 2';
files = dir(folder);
AllDBI = cell(1,1);
AllS = cell(1,1);
AllmedS = cell(1,1);
for m=[24:29]
    if 1
        load(fullfile(files(m).folder,files(m).name)); 
        result = Summary2.result;
        fn = fieldnames(result);
        params = Summary2.params;
        RR = mean(params.experiment.optimalMaps.orig,3)>prctile(reshape(mean(Summary2.params.experiment.optimalMaps.orig,3),[],1),85);
        for j=1
            [~,~,maxind{j}] = retinotopicMapFromIndividualMaps(result.(fn{j}).maps,0,'',93);
            maxind{j} = maxind{j}.*RR;
            clear X;
            for k=1:Summary2.params.experiment.N % iterate the maps of this method
                [row,col] = find(maxind{j}==k);
                X{k,j} = [row,col];
            end
        end
        for j=1
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
            DBI{j} = DBI{j} + mean(TP);
            AllDBI{j} = [AllDBI{j};DBI{j}];
        end % iterate the methods, calculate adjusted silhouette and DBI
    end
end
%% boxplot all
figure;
subplot 121
boxplot([allmeds_loc3;ALLmedS_loc4;ALLmedS_loc8;ALLmedS_9movingbars],makeClustVector([length(allmeds_loc3);length(ALLmedS_loc4);length(ALLmedS_loc8);length(ALLmedS_9movingbars)]),'labels',{'3 moving bars','4 location grid','8 location grid','9 moving bars'},'symbol', '');
ylabel('adjusted Si');box off; set(gca,'XTickLabelRotation',45);
subplot 122
boxplot([alldbi_loc3;AllDBI_loc4;AllDBI_loc8;AllDBI_9movingbars],makeClustVector([length(alldbi_loc3);length(AllDBI_loc4);length(AllDBI_loc8);length(AllDBI_9movingbars)]),'labels',{'3 moving bars','4 location grid','8 location grid','9 moving bars'},'symbol', '');
ylabel('adjusted DBI');box off; set(gca,'XTickLabelRotation',45);
%% main figure 9 - dbi and si of target flankers
folder = 'C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\comparision results 2';
files = dir(folder);
AllDBI = cell(1,1);
AllS = cell(1,1);
AllmedS = cell(1,1);
for m=[24:28]
    if 1
        RR = mean(mapTSCA,3)>prctile(reshape(mean(mapTSCA,3),[],1),93);
        for j=1
%             avgLeftSpatial = mean(mapTSCA(:,:,[1,5,9]),3);
            avgLeftSpatial = mean(mapTSCA(:,:,1:3),3);
%             avgTargetSpatial = mean(mapTSCA(:,:,[2:2:8]),3);
            avgTargetSpatial = mean(mapTSCA(:,:,4:6),3);
%             avgRigthspatial = mean(mapTSCA(:,:,[3,7]),3);
            avgRigthspatial = mean(mapTSCA(:,:,7:9),3);
            [~,~,maxind{j}] = retinotopicMapFromIndividualMaps_target_flanker(cat(3,avgLeftSpatial,avgTargetSpatial,avgRigthspatial),2,'',94);
            maxind{j} = maxind{j}.*RR;
            clear X;
            for k=1:3 % iterate the maps of this method
                [row,col] = find(maxind{j}==k);
                X{k,j} = [row,col];
            end
        end
        for j=1
            [s{j},h] = silhouette(cat(1,X{:,j}),makeClustVector(cellfun(@(x) size(x,1),X(:,j)))');
            NN = histcounts(makeClustVector(cellfun(@(x) size(x,1),X(:,j)))'); % size of each cluster
            expectedClustSize = length(makeClustVector(cellfun(@(x) size(x,1),X(:,j)))')/3; % expected cluster size
            NCM = (3-length(unique(NN))); % number of clusters mismatch
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
    end
end
