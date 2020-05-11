function [status,img,zfit,correlation,fitresult,offsets] = tsca_loc8_2(Z,noises,tlims,settle,gamma,toProject,reduceComp,gaussfiltSTD,toFit,toSave)
% Inputs: parameters for output. to save either empty/zero or contains name
% of experiment (date)
% Outputs: status,signal spatial component, gaussian fit, fitresults,
% offsets
global bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0 signal2_smooth basis
status = 0;
% noise2.time = [mean(cfn{1})]';
% noise2.time = [mean(cfn{2})]';
noise1.time = eye(800)/size(cfn{1},2); % autocorrelation matrix of white noise
if ~isempty(noises)
    for i=1:length(noises)
        noise2(i).time = createToeplitz(noises(i).f0,noises(i).bw,3,[1 0.1 0.1],800);
    end
else
    noise2 = noises;
end

img = cell(8,1); tempimg = cell(8,1);
zfit = cell(8,1);
correlation = zeros(8,1);
fitresult = cell(8,1);
for i=3:length(cfn)
    disp(['loc ' num2str(i-2)]);
    signal.time = [zeros(1,(i-3)*100) [basis(:,1)+basis(:,2)]' zeros(1,(length(cfn)-(i))*100)];
    try
        [ind_mx1]=xplmts(Z, [], 0, [tlims(1) tlims(2)], [settle(1) settle(2)], 0);
    catch
        ind_mx1 = [];
    end
    [projected,components,~,~,output] = tscaFunc(Z,signal,[noise1 noise2],gamma,toProject,reduceComp);
%     tscaAnalyze(output,3,[],1,800);
    estimatedSignal = reshape(components(:,1),size(brn,1),[]);
%     estimatedSignal = MinMaxNorm(estimatedSignal);
    correlation(i-2) = corr(abs(projected(1,:))',signal.time');
%     estimatedSignal = reshape(components(:,1),size(brn0,1),[])+min(min(components(:,1)));
    if ~isempty(ind_mx1)
        estimatedSignal(ind_mx1==0) = 0;
    end
    img{i-2} = reshape(estimatedSignal,size(brn,1),[]);
    if gaussfiltSTD>0
        img{i-2} = imgaussfilt(img{i-2},gaussfiltSTD);
    end
    if toFit
        img2fit = img{i-2}; img2fit = MinMaxNorm(img2fit);img2fit(img2fit<prctile(img2fit(:),85)) = min(img2fit(:));
        [fitresult{i-2}, zfit{i-2}, fiterr1, zerr, resnorm, rr] = fmgaussfit(1:size(brn,2),1:size(brn,1),img2fit);
    end
    estimatedSignal([1:5, end-5:end],:) = min(min(estimatedSignal)); % make edges zero
    estimatedSignal(:,[1:5, end-5:end]) = min(min(estimatedSignal)); % make edges zero
    tempimg{i-2} = reshape(estimatedSignal,size(brn,1),[]);
end
mxc_Ori = zeros(numel(img{3}),length(img)-1);
index2plot = [11,3,15,7,5,13,1,9]; % 8 locs
figure(100);
% suptitle(['gammas: ' num2str(gamma)]); 
% figure(101);
% suptitle(['gammas: ' num2str(gamma)]);
figure(200);
% suptitle('fitted gaussians'); 
figure(300);
% suptitle('waveforms of C.O.M and off C.O.M pixels'); 
t = linspace(0,1,100);
cAxis = [prctile(reshape(cat(3,img{:}),[],1),10) prctile(reshape(cat(3,img{:}),[],1),99)];
for i=1:length(img)
    figure(100);
    subplot(3,5,index2plot(i));
    imf2(img{i}); title(lgn(i+2,:)); 
    caxis(cAxis);
%     caxis([prctile(img{i}(:),10) prctile(img{i}(:),99)]);
    mxc_Ori(:,i) = reshape(img{i},[],1);
    
    figure(200);
    subplot(3,5,index2plot(i));
    imf2(zfit{i}); title(lgn(i+2,:)); colormap(jet);
    
    figure(300);
    subplot(3,5,index2plot(i));
    [Icm] = calcCOM(tempimg{i}-min(tempimg{i}(:)));
%     [~,Icm] = max(reshape(img{i},[],1));
    [Row,Col]= ind2sub(size(brn),Icm);
    [Rows,Cols] = meshgrid(Row-3:Row+3,Col-3:Col+3,length(Row));
    signal = mean(cfn{i+2}(sub2ind(size(brn),Rows,Cols),:));
%     bslnSig = mean(cfn0{1}(sub2ind(size(brn0),Rows,Cols),:));
    % < -- construct baseline signal
    valid_choices = setdiff(1:size(brn,1), Rows);
    nonCOMpxls_Rows = valid_choices(randi(length(valid_choices),length(Rows)));
    valid_choices = setdiff(1:size(brn,2), Cols);
    nonCOMpxls_Cols = valid_choices(randi(length(valid_choices),length(Rows)));
    bslnSig = mean(cfn{i+2}(sub2ind(size(brn),nonCOMpxls_Rows,nonCOMpxls_Cols),:));
    % --->
    [ s_time, s_filteredSignal ] = filter_FIR_and_fix_delay( ones(1,10), 10, t , signal ); % filter to smooth signal
    [ s_timebsln, s_filteredSignalbsln ] = filter_FIR_and_fix_delay( ones(1,10), 10, t , bslnSig ); % filter to smooth signal
    plot(t,signal); hold on; plot(s_time,s_filteredSignal,'linewidth',3); plot(s_timebsln,s_filteredSignalbsln,'linewidth',3);
    title(lgn(i+2,:));xlabel('time [sec]'); ylabel('Amp'); 
%     legend('orig.','smooth','bsln','location','best');
end
% <---- colors map
% [~,I] = max(mxc_Ori,[],2);
% colors = [0 0 0; 0 0 1; 0 1 0; 0 1 1; 1 0 0; 1 0 1; 1 1 0; 1 1 1];
% MAP = colors(I,:);
% figure(400); imagesc(reshape(MAP,size(brn0,1),[],3));
% title('loc 1 - black','loc 2','loc 3','loc 4','loc 5','loc 6','loc 7','loc 8');
% ---->
% <---- calculate offsets
try
    offsets = calcOffsets(fitresult);
catch
    if ~toFit
        disp('user did not specify to fit gaussians');
    else
        disp('could not calculate offsets from the fitted gaussians');
    end
end
% --->
% <---- plot contours
plotContours8(zfit);
% ---->
% <---- save figures and .mat outputs
if toSave
    folder = fullfile('C:\Users\Ori\Desktop\Ori\2nd degree\mtdt',num2str(toSave));
    save(fullfile(folder,'tsca_output'),'img');
%     saveas(gcf,'tscatest_BSLN.jpg');
    save(fullfile(folder,'offests'),'offsets');
end
% ---->
try
    [ds rf cf L]=cmbn84(mxc_Ori, rxc(:,2:end-1), 2); % combine 8 conditions
    status = 1;
catch
    disp('could not combine 8');
end
end

