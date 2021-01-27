function [status,img] = tsca_ff(noises,tlims,settle,gamma,toProject,reduceComp,gaussfiltSTD,toSave)
% Inputs: parameters for output
% Outputs: status,signal spatial component
global bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0 signal2_smooth
status = 0;
noise1.time = eye(100)/size(cfn{1},2); % autocorrelation matrix of white noise
for i=1:length(noises)
    noise2(i).time = createToeplitz(noises(i).f0,noises(i).bw,3,[],100);
end
img = cell(2,1);
signal.time = signal2_smooth;
Means = cellfun(@(x) mean(x,2),cfn0,'uniformoutput',false);
Means2 = zeros(size(cfn0{1},1),1);
for i=1:length(cfn0)
    Means2 = Means2+Means{i};
end
Means2 = Means2./length(Means2);
Z = cfn0{2}-Means2;
try
    [ind_mx1]=xplmts(Z, [], 0, [tlims(1) tlims(2)], [settle(1) settle(2)], 0);
catch
    ind_mx1 = [];
end
[ projected,components] = tscaFunc(Z,signal,[noise1 noise2],gamma,toProject,reduceComp);
estimatedSignal = components(:,1);
estimatedBSLN = components(:,2);
if ~isempty(ind_mx1)
    estimatedSignal(ind_mx1==0) = 0;
end
if ~isempty(gaussfiltSTD)
    img{1} = imgaussfilt(reshape(estimatedSignal,size(brn0,1),[]),gaussfiltSTD);
    img{2} = imgaussfilt(reshape(estimatedBSLN,size(brn0,1),[]),gaussfiltSTD);
end
[img{1}] = fixInversion(img{1});
% [img{2}] = fixInversion(img{2});
figure;
subplot 121; imagesc(img{2}); colormap(jet); title('spatial component - baseline (no flash)');
subplot 122; imagesc(img{1}); colormap(jet); title('spatial component - ff');
figure;
subplot 121; plot(projected(2,:));title('time course - bsln');
subplot 122; plot(projected(1,:));title('time course - ff');
if toSave
    save('tsca_output_new','img');
    saveas(gcf,'tscatest_BSLN.jpg');
end
end

