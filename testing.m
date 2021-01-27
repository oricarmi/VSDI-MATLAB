target_start_frame = [101;401;701];
target_sum = zeros(size(brn));
figure;
for i=1:3
    [Signal1,Signal2,Signal3,beta] = GLM_VSDI(reshape(Z(:,:,target_start_frame(i):target_start_frame(i)+99),[],100),[0.3 3.3 4.5],responseSig');
    subplot(3,1,i);
    imagesc(reshape(beta(end,:),270,327));colormap(gray);
    target_sum = target_sum + reshape(beta(end,:),270,327);
end
figure;imagesc(target_sum)
%%
avgResponse = mean(reshape(c0,[],1000));
gaussEq = 'a*exp(-((x-b)/c)^2)+d';
t = 0:0.01:10-0.01;
fitobj = fit(t',avgResponse',gaussEq);
figure;plot(fitobj,t,avgResponse);
Exponential_Detrend = @(x) fittedmodel.a*exp(fittedmodel.b*x)+fittedmodel.c*exp(fittedmodel.d*x);
exp_t = Exponential_Detrend(t);
detrended_data = avgResponse-exp_t;
figure;plot(t,detrended_data);
exp_t_all = repmat(reshape(exp_t,1,1,[]),270,327);
c_detrended = c0-exp_t_all;
% c_detrended_normalized = (c_detrended-min(c_detrended(:)))/(max(c_detrended(:))-min(c_detrended(:)));
c_detrended_normalized = c_detrended-mean(c_detrended(:,:,1:100),3);
figure
figure;implay(c_detrended_normalized);