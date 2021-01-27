t = 0:0.01:10;
N = length(t);
y = 3+cos(2*pi*2.*t)+normrnd(0,0.5,1,N)...
    +30*responseSig...
    +10*circshift(responseSig,100)...
    +12*circshift(responseSig,200);
y = y';
X = [ones(N,1),cos(2*pi*2.*t)',responseSig',circshift(responseSig,100)',circshift(responseSig,200)'];
figure;plot(t,y);hold on;
beta = (X'*X)^(-1)*X'*y;
plot(t,X*beta);

y_no_noise = y-X(:,1:2)*beta(1:2);
y_no_noise2 = X(:,3:end)*beta(3:end);
figure;
subplot(1,2,1);plot(t,y_no_noise);
subplot(1,2,2);plot(t,y_no_noise2);
error = y-X(:,[1,3:end])*beta([1,3:end]);
figure;
subplot(1,2,1);plot(t,error);
subplot(1,2,2);histogram(error,50);
