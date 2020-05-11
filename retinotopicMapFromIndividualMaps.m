function [retinotopicMap,retinotopicMap2] = retinotopicMapFromIndividualMaps(map,weights,Ttle)
    global brn ump lgn params
    % map should be 3 dimentional, m by m by #signals
    if nargin <3
        Ttle = [];
    end
    if iscell(map)
        map = cat(3,map{:});
    end
    if nargin<2 || ~any(size(weights)==size(map,3))
        [maxmap,maxind] = max(map,[],3); maxmap = maxmap(:); maxind = maxind(:); % get maximum of each pixel, value and index of maximum
    else
        B = bsxfun(@times,map,reshape(weights,1,1,[]));
        [maxmap,maxind] = max(B,[],3); maxmap = maxmap(:); maxind = maxind(:);
    end
    % map22 = colors(maxind,:); map22(maxmap(:)< prctile(maxmap(:),70),:) = repmat([0 0 0],length(find(maxmap(:)< prctile(maxmap(:),70))),1);
    if size(map,1)==40 % simulation
        retinotopicMap = hsv2rgb(maxind/length(unique(maxind)),ones(size(maxind)),(maxmap-min(maxmap))./(max(maxmap)-min(maxmap))); % hue is index of max, saturation is 1 and value is max value (minmax normalized)
        figure;imagesc(reshape(retinotopicMap,size(map,1),size(map,2),3)); %
        figure;imagesc(hsv2rgb([1:length(unique(maxind))]./length(unique(maxind)),ones(1,length(unique(maxind))),ones(1,length(unique(maxind)))));
    else % real data
        retinotopicMap = hsv2rgb(maxind/length(unique(maxind)),ones(size(maxind)),double(maxmap>prctile(maxmap,90)).*(maxmap-min(maxmap))./(max(maxmap)-min(maxmap))); % hue is index of max, saturation is 1 and value is max value (minmax normalized)
        tempbrn = repmat(reshape((brn-min(brn))./max(brn),[],1),1,1,3);
        retinotopicMap2 = squeeze(retinotopicMap); retinotopicMap2(sum(retinotopicMap2,2)==0,:) = tempbrn(sum(retinotopicMap2,2)==0,:); % make brain background instead of black
%             figure;imagesc(ump.*[0:size(brn,2)-1]./1000, ump.*[0:size(brn,1)-1]./1000,reshape(retinotopicMap2,size(map,1),size(map,2),3));title(Ttle); %
        figure;imagesc(reshape(retinotopicMap2,size(map,1),size(map,2),3));title(Ttle); 
        colorz = hsv2rgb([1:length(unique(maxind))]./length(unique(maxind)),ones(1,length(unique(maxind))),ones(1,length(unique(maxind))));
        switch params.experiment.what  % sweep
            case 8
                figure;
                index2plot = [11,3,15,7,5,13,1,9]; % 8 locs
                for i=1:length(colorz)
                    subplot(3,5,index2plot(i));
                    set(gca,'Color',colorz(1,i,:));
                    title(lgn(i+2,:));
                end
            case 9 
                figure;
                index2plot = [5,7,3,8,2,9,1,6,4];
                for i=1:length(colorz)
                    subplot(3,3,index2plot(i));
                    set(gca,'Color',colorz(1,i,:));
                    title(lgn(i+2,:));
                end
            otherwise
                figure;imagesc(colorz); 
        end
    end
end

