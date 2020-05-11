t = 0:1/fs:(1-1/fs);
alphas = [50 200; 10 100; 10 50;50 150]./1000;
randomAlphas = [...
    alphas(1,1)+rand(15000,1)*(alphas(1,2)-alphas(1,1))...
    alphas(2,1)+rand(15000,1)*(alphas(2,2)-alphas(2,1))...
    alphas(3,1)+rand(15000,1)*(alphas(3,2)-alphas(3,1))...
    alphas(4,1)+rand(15000,1)*(alphas(4,2)-alphas(4,1))];
r = zeros(size(randomAlphas,1),length(t));
for i=1:size(r,1) % iterate random signals
    for j=1:size(r,2) % iterate time
        if t(j)<=randomAlphas(i,1) || t(j)>=(randomAlphas(i,1)+randomAlphas(i,2)+randomAlphas(i,3)+randomAlphas(i,4)) % if before latency or after decay, it is zero
            continue
        elseif randomAlphas(i,1)<=t(j) && t(j)<=(randomAlphas(i,1)+randomAlphas(i,2)) % if in rise time
            r(i,j) = 0.5*(1-cos(pi*(t(j)-randomAlphas(i,1))/randomAlphas(i,2)));
        elseif (randomAlphas(i,1)+randomAlphas(i,2))<=t(j) && t(j)<=(randomAlphas(i,1)+randomAlphas(i,2)+randomAlphas(i,3)) % if in plateau time
            r(i,j) = 1;
        else %if (randomAlphas(i,1)+randomAlphas(i,2)+randomAlphas(i,3))<=t(j) && t(j)<=(randomAlphas(i,1)+randomAlphas(i,2)+randomAlphas(i,3)+randomAlphas(i,4)) % if in decaying time
            r(i,j) = 0.5*(1+cos(pi*(t(j)-randomAlphas(i,1)-randomAlphas(i,2)-randomAlphas(i,3))/randomAlphas(i,4)));
        end
    end
end
% figure;plot(r(randomAlphas(:,3)<0.03,:)');
[U,S,V] = svd(r,'econ');
basis = V(:,1:10);
%% non-linear LM
t = [0:0.01:10]';
y1 = 3*t+4 + normrnd(0,1,size(t));
y2 = 3*cos(2*pi*3*t+pi/3) + normrnd(0,1,size(t));
X1 = [ones(length(t),1) t]; 
beta = inv(X1'*X1)*X1'*y1;
yhat1 = X1*beta;
figure;hold on; plot(t,y1); plot(t,yhat1);




    
