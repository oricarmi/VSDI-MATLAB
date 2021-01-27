function [] = saveVideo(cfn,path2save)
figure; 
v = VideoWriter(path2save);
open(v);
zz3D = [];
for i=1:length(cfn)
    zz3D = cat(3,zz3D,rshp(cfn{i}));
end
cAxis = [prctile(zz3D(:),1) prctile(zz3D(:),100)];
for l=1:size(zz3D,3)
    imagesc(zz3D(:,:,l));colormap('gray');
    frame = getframe(gcf);
    writeVideo(v,frame);
end
close(v);
end

