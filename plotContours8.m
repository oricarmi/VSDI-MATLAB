function [fitted2DGContours] = plotContours8(zfits)
% Plot fitted gaussian contours with brain in background
global brn lgn ump
fitted2DGContours = cell(length(zfits)-1,1);
percents = [99.4 99.5; 99.6 99.7; 99.8 99.9; 99.95 100];
resIMG = repmat(reshape((brn-min(brn))./max(brn),[],1),1,1,3);
colors = [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1];
% colors = [1 0 0; 0 1 0; 0 0 1; 1 0 1; 0 1 1; 1 1 0; 1 1 1; 0 0 0];
indContours = zeros(numel(brn),1);
for j= 1:length(zfits)
    zfit = reshape(zfits{j},[],1);
    temp = zeros(size(zfit));
    for i=1:length(percents)
        temp(zfit<prctile(reshape(zfit,[],1),percents(i,2)) & zfit>prctile(reshape(zfit,[],1),percents(i,1))) = 1;
%         temp = zfit<prctile(reshape(zfit,[],1),percents(i,2)) & zfit>prctile(reshape(zfit,[],1),percents(i,1));
    end
    indContours(temp==1) = j;      
%     fitted2DGContours{j-1} = zfit2;  
end
MAP = -1.*ones(length(indContours),3);
MAP(indContours~=0,:) = colors(indContours(indContours~=0),:);
MAP(sum(MAP,2)==-3,:) = resIMG(sum(MAP,2)==-3,:);
figure;
% imagesc(reshape(MAP,size(brn,1),[],3));
imagesc(ump.*[0:size(brn,2)-1]./1000, ump.*[0:size(brn,1)-1]./1000,reshape(MAP,size(brn,1),[],3));

figure;
index2plot = [11,3,15,7,5,13,1,9]; % 8 locs
for i=1:length(colors)
    subplot(3,5,index2plot(i));
    set(gca,'Color',colors(i,:));
    title(lgn(i+2,:));
end
end
