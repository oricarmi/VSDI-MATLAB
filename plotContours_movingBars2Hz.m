function [fitted2DGContours] = plotContours_movingBars2Hz(zfits)
% Plot fitted gaussian contours with brain in background
global brn0 lgn ump brn
fitted2DGContours = cell(length(zfits)-1,1);
percents = [99.4 99.5; 99.6 99.7; 99.8 99.9; 99.95 100];
resIMG = repmat(reshape((brn-min(brn))./max(brn),[],1),1,1,3);
colors = [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1; 0.5 0.5 0.5];
indContours = zeros(numel(brn),1); % to which contour the pixel belongs to (if zero, no contour)
for j= 1:length(zfits)
    zfit = reshape(zfits{j},[],1); 
    temp = zeros(size(zfit)); % temp will be a temporary binary map that indicates if the pixels that belong to this condition's contours
    for i=1:length(percents)
        temp(zfit<prctile(reshape(zfit,[],1),percents(i,2)) & zfit>prctile(reshape(zfit,[],1),percents(i,1))) = 1;
%         temp = zfit<prctile(reshape(zfit,[],1),percents(i,2)) & zfit>prctile(reshape(zfit,[],1),percents(i,1));
    end
    indContours(temp==1) = j;      
%     fitted2DGContours{j-1} = zfit2;  
end
MAP = -1.*ones(length(indContours),3);
MAP(indContours~=0,:) = colors(indContours(indContours~=0),:); % if there is a contour, put in that line the corresponding rgb color
MAP(sum(MAP,2)==-3,:) = resIMG(sum(MAP,2)==-3,:); % where there isn't contour, put res image pixel
figure;
imagesc(ump.*[0:size(brn,2)-1]./1000, ump.*[0:size(brn,1)-1]./1000,reshape(MAP,size(brn,1),[],3));

figure;
for i=1:length(colors)
    subplot(3,3,i);
    set(gca,'Color',colors(i,:));
    title(lgn(i+3,:));
end
end
