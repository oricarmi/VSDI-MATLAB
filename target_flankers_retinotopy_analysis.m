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
imf2(rshp(tempbrn)); hold on; imf2(avgLeftSpatial,0.85);
title('target');
subplot(1,3,2)
% avgTargetSpatial = mean(mapTSCA(:,:,[2:8]),3);
avgTargetSpatial = mean(mapTSCA(:,:,4:6),3);
imf2(rshp(tempbrn)); hold on; imf2(avgTargetSpatial,0.85);
title('flankers');
subplot(1,3,3)
% avgRigthspatial = mean(mapTSCA(:,:,[3,7]),3);
avgRigthspatial = mean(mapTSCA(:,:,7:9),3);
imf2(rshp(tempbrn)); hold on; imf2(avgRigthspatial,0.85);
title('combined');
figure;
[~,r] = retinotopicMapFromIndividualMaps_target_flanker(cat(3,avgLeftSpatial,avgTargetSpatial,avgRigthspatial),1,'',96);
%%
fn = fieldnames(result);
retmapsAll = figure("name","retmapsAll",'Position', [100 0 1400 500]);
% indmapsAll = figure("name","indmapsAll");
% Titles = {'(1)','(2)','(3)','(4)','(5)','(6)','(7)','(8)'};
Titles ={'TSCA','Tmax','AOF','Corr','GLM','MPT','TSCA and GLM'};
flag=0;
for i=1:6
    thisMap = result.(fn{i}).maps;
    avgLeftSpatial = mean(thisMap(:,:,[1,5,9]),3);
    avgTargetSpatial = mean(thisMap(:,:,[2:8]),3);
    avgRigthspatial = mean(thisMap(:,:,[3,7]),3);
    [~,r] = retinotopicMapFromIndividualMaps_target_flanker(cat(3,avgLeftSpatial,avgTargetSpatial,avgRigthspatial),0,fn{i},94);
    figure(retmapsAll); subplot(2,4,i)
    imf2(r);xlabel('');set(gca,'xticklabel',[]);set(gca,'xtick',[]);ylabel('');set(gca,'yticklabel',[]);set(gca,'ytick',[]);
    title(Titles{i});
end

% avgLeftSpatial = mean(mapTSCA(:,:,[1,5,9]),3);
% avgTargetSpatial = mean(mapTSCA(:,:,[2:8]),3);
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
X = cell(params.experiment.N,1); % preallocate memory
disp('select regoin of interest and double click inside when finished');
figure;RR = roipoly(rshp(retinotopicMap));
maxind = double(RR).*maxind;
for k=1:p3 % iterate the maps 
    [row,col] = find(maxind==k); % get pixels of this map (condition)
    X = [row,col];
end
figure;[s,h] = silhouette(cat(1,X{:}),makeClustVector(cellfun(@(x) size(x,1),X))');
NN = histcounts(makeClustVector(cellfun(@(x) size(x,1),X))');
expectedClustSize = length(makeClustVector(cellfun(@(x) size(x,1),X))')/3;
NNN = abs(NN - expectedClustSize)/expectedClustSize; 
NNN = NNN + (8-length(unique(NN)));
NNNN = repelem(NNN,NN);
s = s - NNNN';
s = 2*(s-min(s))/(max(s)-min(s))-1;
[~,~,~,~,~,DBI] = ClusterSimilarity(X);