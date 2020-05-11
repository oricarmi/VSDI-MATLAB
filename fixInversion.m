function [imgFixed] = fixInversion(img,projected,signal)
% Check if inputted image is inverted 
%     crossRows = 110:160;
%     crossCols = 150:180;
%     if sum((projected'-signal.time).^2)>sum(((-1).*projected'-signal.time).^2) % inverted time approximation, so invert
%         img = img.*(-1);
%     end
    imgFixed = sign(projected(1,:)*signal').*img; % multipy by sign of inner product to fix inversion
%     img = sign(corr(projected(1,:)',signal.time)).*img; % multipy by sign of inner product to fix inversion
%     if length(find(img>0))>length(find(img<0)) % && abs(max(reshape(img,[],1)))<abs(min(reshape(img,[],1))) 
%         img = -1.*img;
%     end
%     if mean(mean(img(crossRows,crossCols)))<mean(mean(img(setdiff(1:size(brn,1),crossRows),setdiff(1:size(brn,2),crossCols))))
%         img = -1.*img;
%     end
%     if sum(fiterr1)>sum(fiterr2)
%         img{i-1} = -1.*img{i-1};
%         zfit{i-1} = zfit2;
%     end
%% Version 1
%     imgFixed = img;
%     wantedMap1 = zeros(size(img,1),size(img,2));
%     wantedMap2 = zeros(size(img,1),size(img,2));
%     range = abs(max(reshape(img,[],1))-min(reshape(img,[],1)));
%     [X,Y] = meshgrid(1:size(wantedMap1,2),1:size(wantedMap1,1));
% %     cX = size(wantedMap,2)/2; cY = size(wantedMap,1)/2; 
%     r = 15;
%     [~,cX] = max(max(img)); [~,cY] = max(max(img,[],2));
%     wantedMap1((X-cX).^2+(Y-cY).^2<r.^2) = range;
%     [~,cX] = max(max(-1.*img)); [~,cY] = max(max(-1.*img,[],2));
%     wantedMap2((X-cX).^2+(Y-cY).^2<r.^2) = range;
%     C1 = xcorr2(img,wantedMap1); 
%     C2 = xcorr2(-1.*img,wantedMap2);
%     if max(max(C2))>max(max(C1))
%         imgFixed = -1.*img;
%     end
%% Version 2
%     imgFixed = img;
%     wantedMap = zeros(size(img));
%     range = abs(max(reshape(img,[],1))-min(reshape(img,[],1)));
%     [X,Y] = meshgrid(1:size(wantedMap,2),1:size(wantedMap,1));
% %     cX = size(wantedMap,2)/2; cY = size(wantedMap,1)/2; 
%     r = 15;
%     [~,cX] = max(max(img)); [~,cY] = max(max(img,[],2));
%     wantedMap((X-cX).^2+(Y-cY).^2<r.^2) = range;
%     C1 = xcorr2(img,wantedMap); 
%     C2 = xcorr2(-1.*img,wantedMap);
%     if max(max(C2))>max(max(C1))
%         imgFixed = -1.*img;
%     end
end

