%% imshow the general responses obtained by tsca&glm
ZZ3d = rshp(ZZ);
figure;
tempbrn = repmat(reshape((brn-min(brn))./(max(brn)-min(brn)),[],1),1,1,3);
for i=1:9
    subplot(3,3,i)
    imf2(rshp(tempbrn)); hold on; imf2(mapTSCA(:,:,i),0.85);
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
% avgLeftSpatial = mean(mapTSCA(:,:,[1,5,9]),3);
avgLeftSpatial = mean(mapTSCA(:,:,1:3),3);
imf2(rshp(tempbrn)); hold on; imf2(avgLeftSpatial,0.55);
title('target');
subplot(1,3,2)
% avgTargetSpatial = mean(mapTSCA(:,:,[2:2:8]),3);
avgTargetSpatial = mean(mapTSCA(:,:,4:6),3);
imf2(rshp(tempbrn)); hold on; imf2(avgTargetSpatial,0.75);
title('flankers');
subplot(1,3,3)
% avgRigthspatial = mean(mapTSCA(:,:,[3,7]),3);
avgRigthspatial = mean(mapTSCA(:,:,7:9),3);
imf2(rshp(tempbrn)); hold on; imf2(avgRigthspatial,0.55);
title('combined');
figure;
[~,~,maxind] = retinotopicMapFromIndividualMaps_target_flanker(cat(3,avgLeftSpatial,avgTargetSpatial,avgRigthspatial),2,'',80);
%%
fn = fieldnames(result);fn = fn([7,1:6]);
retmapsAll = figure("name","retmapsAll",'Position', [100 0 1400 500]);
% indmapsAll = figure("name","indmapsAll");
% Titles = {'(1)','(2)','(3)','(4)','(5)','(6)','(7)','(8)'};
Titles ={'TSCA and GLM','TSCA','Tmax','AOF','Corr','GLM','MPT'};
flag=0;
for i=1:7
    thisMap = result.(fn{i}).maps;
    figure(retmapsAll); subplot(2,4,i);
%     avgLeftSpatial = mean(thisMap(:,:,[1,5,9]),3);
%     avgTargetSpatial = mean(thisMap(:,:,[2:2:8]),3);
%     avgRigthspatial = mean(thisMap(:,:,[3,7]),3);
    avgLeftSpatial = mean(thisMap(:,:,[1:3]),3);
    avgTargetSpatial = mean(thisMap(:,:,[4:6]),3);
    avgRigthspatial = mean(thisMap(:,:,[7:9]),3);
    [~,r] = retinotopicMapFromIndividualMaps_target_flanker(cat(3,avgLeftSpatial,avgTargetSpatial,avgRigthspatial),0,fn{i},94);
    imf2(r);
    xlabel('');set(gca,'xticklabel',[]);set(gca,'xtick',[]);ylabel('');set(gca,'yticklabel',[]);set(gca,'ytick',[]);
    title(Titles{i});
end

% avgLeftSpatial = mean(mapTSCA(:,:,[1,5,9]),3);
% avgTargetSpatial = mean(mapTSCA(:,:,[2:2:8]),3);
% avgRigthspatial = mean(mapTSCA(:,:,[3,7]),3);
avgLeftSpatial = mean(mapTSCA(:,:,1:3),3);
avgTargetSpatial = mean(mapTSCA(:,:,4:6),3);
avgRigthspatial = mean(mapTSCA(:,:,7:9),3);
figure(retmapsAll); subplot(2,4,7);
[~,r] = retinotopicMapFromIndividualMaps_target_flanker(cat(3,avgLeftSpatial,avgTargetSpatial,avgRigthspatial),0,'',94);
imf2(r);xlabel('');set(gca,'xticklabel',[]);set(gca,'xtick',[]);ylabel('');set(gca,'yticklabel',[]);set(gca,'ytick',[]);
title(Titles{end}); 
%%
[~,retinotopicMap,maxind] = retinotopicMapFromIndividualMaps_target_flanker(cat(3,avgLeftSpatial,avgTargetSpatial,avgRigthspatial),0,'',96);
X = cell(3,1); % preallocate memory
disp('select regoin of interest and double click inside when finished');
figure;RR = roipoly(rshp(retinotopicMap));
maxind = double(RR).*maxind;
for k=1:3 % iterate the maps 
    [row,col] = find(maxind==k); % get pixels of this map (condition)
    X = [row,col];
end
figure;[s,h] = silhouette(cat(1,X{:}),makeClustVector(cellfun(@(x) size(x,1),X))');
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