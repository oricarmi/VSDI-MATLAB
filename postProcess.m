function [map2] = postProcess(map)
% Post process the maps generated
global params fs 
    map2 = rshp(map);
    if all(params.post.medFiltSize) && params.post.gaussfltSTD
        for i=1:size(map2,3)
            map2(:,:,i) = medfilt2(map2(:,:,i),[params.post.medFiltSize(1) params.post.medFiltSize(2)]);
            map2(:,:,i) = imgaussfilt(map2(:,:,i),params.post.gaussfltSTD);
            switch params.post.normalization
                case 'z'
                    tmp = reshape(map2(:,:,i),[],1);
                    map2(:,:,i) = (map2(:,:,i)-mean(tmp))./std(tmp);
                case 'minmax'
                    map2(:,:,i) = MinMaxNorm(map2(:,:,i));
                otherwise
                    ;
            end
        end
    elseif all(params.post.medFiltSize) % only meidan filter
        for i=1:size(map2,3)
            map2(:,:,i) = medfilt2(map2(:,:,i),[params.post.medFiltSize(1) params.post.medFiltSize(2)]);
            switch params.post.normalization
                case 'z'
                    tmp = reshape(map2(:,:,i),[],1);
                    map2(:,:,i) = (map2(:,:,i)-mean(tmp))./std(tmp);
                case 'minmax'
                    map2(:,:,i) = MinMaxNorm(map2(:,:,i));
                otherwise
                    ;
            end

        end
    elseif params.post.gaussfltSTD % only gaussfilt
        for i=1:size(map2,3)
            map2(:,:,i) = imgaussfilt(map2(:,:,i),params.post.gaussfltSTD);
            switch params.post.normalization
                case 'z'
                    tmp = reshape(map2(:,:,i),[],1);
                    map2(:,:,i) = (map2(:,:,i)-mean(tmp))./std(tmp);
                case 'minmax'
                    map2(:,:,i) = MinMaxNorm(map2(:,:,i));
                otherwise
                    ;
            end
        end
    end

end

