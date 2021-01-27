function [HD] = mvn_Hellinger_distance(mu1,mu2,sigma1,sigma2)
% calculate multi-variate normal distribution's hellinger distance
    if isrow(mu1)
        mu1 = mu1';
    end
    if isrow(mu2)
        mu2 = mu2';
    end
    HD = 1 - ( (det(sigma1)^0.25*det(sigma2)^0.25) / (det((sigma2+sigma1)/2)^0.5) ) * exp( -0.125 * (mu1-mu2)' * ((sigma2+sigma1)/2)^(-1) * (mu1-mu2) );
end

