function [Icm,R,C] = calcCOM(img)
% Calculate center of mass of an image
    global brn
    if isvector(img)
        img = reshape(img,size(brn,1),[]);
    end
    tot_mass = sum(img(:));
    [ii,jj] = ndgrid(1:size(img,1),1:size(img,2));
    R = round(sum(ii(:).*img(:))/tot_mass);
    C = round(sum(jj(:).*img(:))/tot_mass);
    Icm = sub2ind(size(img),R,C);
end

