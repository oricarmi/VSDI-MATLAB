function [Signal1,Signal2,Signal3,beta,Noise] = GLM_VSDI(Z,noiseFreqs,basis)
% Fit a generalizd linear model to the data according to the noises
% specified in noises 
    global fs brn
    T = size(Z,2); % total time
    Z = Z';
    t = [linspace(0,(T-1)./fs,T)]'; % time vector
    X1 = ones(T,1); % DC
    Xn = zeros(T,length(noiseFreqs)*2);
    k = 1;
    for i=1:length(noiseFreqs)
        Xn(:,k:k+1) = [sin(2*pi*noiseFreqs(i).*t) cos(2*pi*noiseFreqs(i).*t)];
        k = k+2;
    end
    Xtot = [X1 Xn basis];
    beta = inv(Xtot'*Xtot)*Xtot'*Z;
    All = Xtot*beta;
    Residuals = Z - All;
    Signal1 = [1/beta(1)*sign(beta(1))*(Z - Xn*beta(2:(size(Xn,2)+1),:))]'; % all signal minus regressors of noise
    Signal2 = [(basis*beta(end-(size(basis,2)-1):end,:)+Residuals)]'; % signal regressor plus residuals
    Signal3 = [basis*beta(end-(size(basis,2)-1):end,:)]'; % signal regressor
    Noise = [Xn*beta(2:(size(Xn,2)+1),:)]'; % noise regressors  
end

