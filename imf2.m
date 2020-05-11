function im_out=imf2(im)
%% show scaled image
global sz rot ump

if size(im,2)==1
    im_out=rot90(reshape(im, sz(1), sz(2)),rot);
else
    im_out=rot90(im ,rot);
end
imagesc(ump.*[0:sz(2)-1]./1000, ump.*[0:sz(1)-1]./1000, im_out);
