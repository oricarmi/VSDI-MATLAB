function [status,img,zfit,fitresult,offests] = tsca_loc9(noises,tlims,settle,gamma,toProject,reduceComp,gaussfiltSTD,toFit,toSave)
% Inputs: parameters for output
% Outputs: status,signal spatial component
global bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0 signal2_smooth
status = 0;
% noise2.time = [mean(cfn{1})]';
% noise2.time = [mean(cfn{2})]';
noise1.time = eye(100)/size(cfn{1},2); % autocorrelation matrix of white noise
for i=1:length(noises)
    noise2(i).time = createToeplitz(noises(i).f0,noises(i).bw,3);
end
Means = cellfun(@(x) mean(x,2),cfn0,'uniformoutput',false);
Means2 = zeros(size(cfn0{1},1),1);
for i=1:length(cfn0)
    Means2 = Means2+Means{i};
end
Means2 = Means2./length(Means2);
img = cell(9,1);
zfit = cell(9,1);
fitresult = cell(9,1);
for i=2:length(cfn0)
    disp(['loc ' num2str(i-1)]);
    signal.time = signal2_smooth;
    Z = cfn0{i}-Means2;
    try
        [ind_mx1]=xplmts(Z, [], 0, [tlims(1) tlims(2)], [settle(1) settle(2)], 0);
    catch
        ind_mx1 = [];
    end
    [ projected,components] = tscaFunc(Z,signal,[noise1 noise2],gamma,toProject,reduceComp);
    estimatedSignal = reshape((components(:,1)-min(components(:,1)))/(max(components(:,1))-min(components(:,1))),size(brn0,1),[]);
%     estimatedSignal = reshape(components(:,1),size(brn0,1),[]);
    estimatedSignal([1:2, end-2:end],:) = min(min(estimatedSignal)); % make edges zero
    estimatedSignal(:,[1:2, end-2:end]) = min(min(estimatedSignal)); % make edges zero
    if ~isempty(ind_mx1)
        estimatedSignal(ind_mx1==0) = min(min(estimatedSignal));
    end
    if ~isempty(gaussfiltSTD)
        img{i-1} = imgaussfilt(estimatedSignal,gaussfiltSTD);
    end
    img{i-1} = fixInversion(img{i-1},projected,signal2_smooth);
    if toFit
        [fitresult{i-1}, zfit{i-1}, fiterr1, zerr, resnorm, rr] = fmgaussfit(1:size(brn0,2),1:size(brn0,1),img{i-1});
    end
end 
mxc_Ori = zeros(numel(img{3}),length(img)-1);
index2plot = [5,7,3,8,2,9,1,6,4];
figure(100);
suptitle(['gammas: ' num2str(gamma)]); 
% figure(101);
% suptitle(['gammas: ' num2str(gamma)]);
figure(200);
suptitle('fitted gaussians'); 
figure(300);
suptitle('waveforms of C.O.M and off C.O.M pixels'); 
t = linspace(0,1,100);
for i=1:length(img)
    figure(100);
    subplot(3,3,index2plot(i));
    imagesc(img{i}); title(lgn(i+1,:)); colormap(jet);caxis([prctile(img{i}(:),5) prctile(img{i}(:),99)]);
    mxc_Ori(:,i) = reshape(img{i},[],1);
    
    figure(200);
    subplot(3,3,index2plot(i));
    imagesc(zfit{i}); title(lgn(i+1,:)); colormap(jet);
    
    figure(300);
    subplot(3,3,index2plot(i));
%     [Icm] = calcCOM(img{i});
    [~,Icm] = max(reshape(img{i},[],1));
    [Row,Col]= ind2sub(size(brn0),Icm);
    [Rows,Cols] = meshgrid(Row-3:Row+3,Col-3:Col+3,length(Row));
    signal = mean(cfn0{i+1}(sub2ind(size(brn0),Rows,Cols),:));
%     bslnSig = mean(cfn0{1}(sub2ind(size(brn0),Rows,Cols),:));
    % < -- construct baseline signal
    valid_choices = setdiff(1:size(brn0,1), Rows);
    nonCOMpxls_Rows = valid_choices(randi(length(valid_choices),length(Rows)));
    valid_choices = setdiff(1:size(brn0,2), Cols);
    nonCOMpxls_Cols = valid_choices(randi(length(valid_choices),length(Rows)));
    bslnSig = mean(cfn0{i+1}(sub2ind(size(brn0),nonCOMpxls_Rows,nonCOMpxls_Cols),:));
    % --->
    % < --- filter & plot the signals
    [ s_time, s_filteredSignal ] = filter_FIR_and_fix_delay( ones(1,10), 10, t , signal ); % filter to smooth signal
    [ s_timebsln, s_filteredSignalbsln ] = filter_FIR_and_fix_delay( ones(1,10), 10, t , bslnSig ); % filter to smooth signal
    plot(t,signal); hold on; plot(s_time,s_filteredSignal,'linewidth',3); plot(s_timebsln,s_filteredSignalbsln,'linewidth',3);
    title(lgn(i+1,:));xlabel('time [sec]'); ylabel('Amp'); 
%     legend('orig.','smooth','bsln','location','best');
    % --->
    
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
plotContours9(zfit);
% ---->
% <---- save figures and .mat outputs
if toSave
    save('tsca_output_new','img');
    saveas(gcf,'tscatest_BSLN.jpg');
end
try
    pc=0; [ry cx mx]=cmbn9(mxc_Ori); 
    status = 1;
catch
    disp('could not combine 9');
end
end

