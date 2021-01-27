function [Z] = showVideo(cfn)
% input: nothing because all the needed variables are global
% output: Z - frames one after the other
global bsln fs sz ump rot fgn brn frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms
Z = [];
for i=1:length(cfn)
    Z = cat(3,Z,rshp(cfn{i}));
end
figure;implay(Z);

