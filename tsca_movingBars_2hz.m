function [status,img,correlation,zfit,fitresult,offsets] = tsca_movingBars_2hz(Z,noises,tlims,settle,gamma,toProject,reduceComp,gaussfiltSTD,toFit,toSave)
% Inputs: parameters for output
% Outputs: status,signal spatial component
global bsln fs sz ump rot fgn brn brn0 frq cmap lgn scl lgn00 fnm c_f vms plt_on pt pc vc xs prms cfn cfn0 signal2_smooth
status = 0;
% noise2.time = [mean(cfn{1})]';
% noise2.time = [mean(cfn{2})]';
noise1.time = eye(size(Z,2))/size(cfn{1},2); % autocorrelation matrix of white noise
for i=1:length(noises)
    noise2(i).time = createToeplitz(noises(i).f0,noises(i).bw,3,[1 0.1 0.1],size(Z,2));
end

img = cell(9,1);
zfit = cell(9,1);
fitresult = cell(9,1);
correlation = zeros(1,9);
for i=4:length(cfn)
    disp(['sweep ' num2str(i-3)]);
    signal.time = [zeros(1,(i-4)*50) signal2_smooth zeros(1,(length(cfn)-(i))*50)];
    try
        [ind_mx1]=xplmts(Z, [], 0, [tlims(1) tlims(2)], [settle(1) settle(2)], 0);
    catch
        ind_mx1 = [];
    end
    [projected,components,~,~,output] = tscaFunc(Z,signal,[noise1 noise2],gamma,toProject,reduceComp);
%     tscaAnalyze(output,3,[],1,size(Z,2));
%     [maxCorrelation,indmaxCorrelation] = max(corr(abs(projected'),signal.time')); % get the index of component with highest correlation to original signal
    correlation(i-2) = corr(abs(projected(1,:))',signal.time');
    estimatedSignal = reshape(components(:,1),size(brn,1),[]);
%     estimatedSignal = MinMaxNorm(estimatedSignal);
    estimatedSignal([1:2, end-2:end],:) = min(min(estimatedSignal)); % make edges zero
    estimatedSignal(:,[1:2, end-2:end]) = min(min(estimatedSignal)); % make edges zero
    if ~isempty(ind_mx1)
        estimatedSignal(ind_mx1==0) = min(min(estimatedSignal));
    end   
    if gaussfiltSTD>0
        img{i-3} = imgaussfilt(estimatedSignal,gaussfiltSTD);
    else
        img{i-3} = estimatedSignal;
    end
%     img{i-1} = fixInversion(img{i-1},projected,signal2_smooth); 
    if toFit
        img2fit = img{i-3}; img2fit(img2fit<prctile(img2fit(:),85)) = min(img2fit(:));
        [fitresult{i-3}, zfit{i-3}, fiterr1, zerr, resnorm, rr] = fmgaussfit(1:size(brn,2),1:size(brn,1),img2fit);
    end
end 
mxc_Ori = zeros(numel(img{3}),length(img)-1);
figure(100);
% suptitle(['gammas: ' num2str(gamma)]); 
% figure(101);
% suptitle(['gammas: ' num2str(gamma)]);
figure(200);
% suptitle('fitted gaussians'); 
figure(300);
% suptitle('waveforms of C.O.M and off C.O.M pixels'); 
t = linspace(0,0.5,50);
cAxis = [prctile(reshape(cat(3,img{:}),[],1),10) prctile(reshape(cat(3,img{:}),[],1),99)];
for i=1:length(img)
    % < ---- plot estimated maps by tsca
    figure(100);
    subplot(3,3,i);
    imf2(img{i}); title(lgn(i+3,:)); colormap(jet);
%     caxis([prctile(img{i}(:),10) prctile(img{i}(:),99)]);
    caxis(cAxis);
    mxc_Ori(:,i) = reshape(img{i},[],1);
    % --->
    % <--- plot fitted gaussians
    figure(200);
    subplot(3,3,i);
    imf2(zfit{i}); title(lgn(i+3,:)); colormap(jet);
    % --->
    % <---- plot time course of on and off C.O.M pixels (response and no response)
    figure(300);
    subplot(3,3,i);
%     [Icm] = calcCOM(img{i});
    tempIMG = img{i};
    tempIMG([1:2, end-2:end],:) = min(min(tempIMG)); % make edges zero
    tempIMG(:,[1:2, end-2:end]) = min(min(tempIMG)); % make edges zero
    [~,Icm] = max(reshape(tempIMG,[],1));
    [Row,Col]= ind2sub(size(brn),Icm);
    [Rows,Cols] = meshgrid(Row-3:Row+3,Col-3:Col+3,length(Row));
    signal = mean(cfn{i+3}(sub2ind(size(brn),Rows,Cols),:));
%     bslnSig = mean(cfn0{1}(sub2ind(size(brn0),Rows,Cols),:));
    % < -- construct baseline signal (from off C.O.M pixels)
    valid_choices = setdiff(1:size(brn,1), Rows);
    nonCOMpxls_Rows = valid_choices(randi(length(valid_choices),length(Rows)));
    valid_choices = setdiff(1:size(brn,2), Cols);
    nonCOMpxls_Cols = valid_choices(randi(length(valid_choices),length(Rows)));
    bslnSig = mean(cfn{i+3}(sub2ind(size(brn),nonCOMpxls_Rows,nonCOMpxls_Cols),:));
    % --->
    % < --- filter & plot the signals
    [ s_time, s_filteredSignal ] = filter_FIR_and_fix_delay( ones(1,10), 10, t , signal ); % filter to smooth signal
    [ s_timebsln, s_filteredSignalbsln ] = filter_FIR_and_fix_delay( ones(1,10), 10, t , bslnSig ); % filter to smooth signal
    plot(t,signal); hold on; plot(s_time,s_filteredSignal,'linewidth',3); plot(s_timebsln,s_filteredSignalbsln,'linewidth',3);
    title(lgn(i+3,:));xlabel('time [sec]'); ylabel('Amp'); 
%     legend('orig.','smooth','bsln','location','best');
    % --->
    % ----->
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
    offsets = [];
    if ~toFit
        disp('user did not specify to fit gaussians');
    else
        disp('could not calculate offsets from the fitted gaussians');
    end
end
% --->
% <---- plot contours
plotContours_movingBars2Hz(zfit);
% ---->
% <---- save figures and .mat outputs
if toSave
    save('tsca_output_new','img');
    saveas(gcf,'tscatest_BSLN.jpg');
end
try
    pc=0; [ry cx mx]=cmbn9_movingBars(mxc_Ori); 
    status = 1;
catch
    disp('could not combine 9');
end
end

