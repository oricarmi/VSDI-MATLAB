function plotMaps(map,Ttle,CA)
% plot maps (3d matrix, 3rd dimension is differnet condition)
global params lgn
    if nargin<2
        CA = 0;
        Ttle = 'unknown';
    elseif nargin<3
        CA = 0;
    end
    if iscell(map)
        map = cat(3,map{:});
    end 
    figure("name",sprintf('indMaps %s', Ttle)); 
    suptitle([Ttle]);% ", caxis: " num2str(CA)]);
    cAxis = [prctile(map(:),1) prctile(map(:),100)];
    switch params.experiment.what
        case 8
            index2plot = [11,3,15,7,5,13,1,9]; % 8 locs
            for i=1:8
                subplot(3,5,index2plot(i));
                imagesc(map(:,:,i));title(lgn(i+2,:));
                if CA
                    caxis(cAxis);
                end
            end
        case 9
            index2plot = [5,7,3,8,2,9,1,6,4];
            for i=1:9
                subplot(3,3,index2plot(i));
                imagesc(map(:,:,i));title(lgn(i+1,:));
                if CA
                    caxis(cAxis);
                end
            end
        case 4
            index2plot = [3,2,1,4];
            for i=1:4
                subplot(2,2,index2plot(i));
                imagesc(map(:,:,i));title(lgn(i+2,:));
                if CA
                    caxis(cAxis);
                end
            end
        case 52
            index2plot = [7,3,1,9,5];
            for i=1:5
                subplot(3,3,index2plot(i));
                imagesc(map(:,:,i));title(lgn(i+2,:));
                if CA
                    caxis(cAxis);
                end
            end
        otherwise % moving bars 2[Hz]
            for i=1:9
                subplot(3,3,i)
                imagesc(map(:,:,i));title(lgn(i+3,:));
                if CA
                    caxis(cAxis);
                end
            end
    end

end

