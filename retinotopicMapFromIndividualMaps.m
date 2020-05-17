function [retinotopicMap,retinotopicMap2] = retinotopicMapFromIndividualMaps(maps,weights,Ttle,pTile,manualSelection)
    global brn ump lgn params
    % map should be 3 dimentional, m by m by #signals
    switch nargin
        case 4
            manualSelection = 0;
        case 3
            pTile = 90;
            manualSelection = 0;
        case 2
            Ttle = [];
            pTile = 90;
            manualSelection = 0;  
    end
    if iscell(maps)
        maps = cat(3,maps{:});
    end
    if manualSelection
        tempmaps = maps; clear maps; 
        for i=1:size(tempmaps,3) % min max norm the maps so that maximum is the same height (scaling)
            maps(:,:,i) = tempmaps(manualSelection(2):manualSelection(2)+manualSelection(4),manualSelection(1):manualSelection(1)+manualSelection(3),i);
            maps(:,:,i) = MinMaxNorm(maps(:,:,i));   
        end
        brnn = brn(manualSelection(2):manualSelection(2)+manualSelection(4),manualSelection(1):manualSelection(1)+manualSelection(3));
        tempbrn = repmat(reshape((brnn-min(brnn))./(max(brnn)-min(brnn)),[],1),1,1,3);
    else
        for i=1:size(maps,3) % min max norm the maps so that maximum is the same height (scaling)
            maps(:,:,i) = MinMaxNorm(maps(:,:,i));   
        end
        tempbrn = repmat(reshape((brn-min(brn))./(max(brn)-min(brn)),[],1),1,1,3);
    end
    if nargin<2 || ~any(size(weights)==size(maps,3)) % if there arent enough weights
        [maxmap,maxind] = max(maps,[],3); maxmap = maxmap(:); maxind = maxind(:); % get maximum of each pixel (value and index)
    else
        B = bsxfun(@times,maps,reshape(weights,1,1,[])); % multiply each map by its weight
        [maxmap,maxind] = max(B,[],3); maxmap = maxmap(:); maxind = maxind(:); % get maximum of each pixel (value and index) after weighting
    end
    if size(maps,1)==40 % simulation
        retinotopicMap = hsv2rgb(maxind/length(unique(maxind)),ones(size(maxind)),(maxmap-min(maxmap))./(max(maxmap)-min(maxmap))); % hue is index of max, saturation is 1 and value is max value (minmax normalized)
        figure;imagesc(reshape(retinotopicMap,size(maps,1),size(maps,2),3)); %
        figure;imagesc(hsv2rgb([1:length(unique(maxind))]./length(unique(maxind)),ones(1,length(unique(maxind))),ones(1,length(unique(maxind)))));
    else % real data
        retinotopicMap = hsv2rgb(maxind/length(unique(maxind)),ones(size(maxind)),double(maxmap>prctile(maxmap,pTile)).*maxmap); % hue is index of max, saturation is 1 and value is max value (minmax normalized)
        retinotopicMap2 = squeeze(retinotopicMap); retinotopicMap2(sum(retinotopicMap2,2)==0,:) = tempbrn(sum(retinotopicMap2,2)==0,:); % make brain background instead of black
        figure;imagesc(ump.*[0:size(brn,2)-1]./1000, ump.*[0:size(brn,1)-1]./1000,reshape(retinotopicMap2,size(maps,1),size(maps,2),3));title(Ttle); %
%         figure;imagesc(reshape(retinotopicMap2,size(map,1),size(map,2),3));title(Ttle); 
        colorz = hsv2rgb([1:length(unique(maxind))]./length(unique(maxind)),ones(1,length(unique(maxind))),ones(1,length(unique(maxind))));
        switch params.experiment.what  
            case 8 % loc 8
                figure;
                index2plot = [11,3,15,7,5,13,1,9]; % 8 locs
                for i=1:length(colorz)
                    subplot(3,5,index2plot(i));
                    set(gca,'Color',colorz(1,i,:));
                    title(lgn(i+2,:));
                end
            case 9  % loc 9
                figure;
                index2plot = [5,7,3,8,2,9,1,6,4];
                for i=1:length(colorz)
                    subplot(3,3,index2plot(i));
                    set(gca,'Color',colorz(1,i,:));
                    title(lgn(i+1,:));
                end
            case 4 % loc 4
                figure;
                index2plot = [3,2,1,4];
                for i=1:length(colorz)
                    subplot(2,2,index2plot(i));
                    set(gca,'Color',colorz(1,i,:));
                    title(lgn(i+2,:));
                end
            otherwise % sweep
                figure;imagesc(colorz); 
        end
    end
end

