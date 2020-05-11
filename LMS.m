function [ Clean_Signal,W ] = LMS( Raw_Signal , Noise_ref, Options )
% A function that cleans the signal according to the noise reference
% sequences with LMS algorithm
    N = length(Raw_Signal); % get the length of the signal;
    % < ---- make input row signals, as the function works for row signal inputs       
    if ~iscolumn(Raw_Signal)
        Raw_Signal = Raw_Signal';
    end
    if size(Noise_ref,1) ~= N
        Noise_ref = Noise_ref';
    end
    % --->
    R = Noise_ref'*Noise_ref/N; % get the expected value (mean) of each element in R (size is M x M)
    M = size(R,1); % number of reference signals
    % <--- define default W0 and mu if not sent as argument to function 
    if nargin<3
        mu = 0.5/trace(R);
        W0 = zeros(M,1);
    else
        if Options.mu == 0 || ~isfield(Options,'mu')
            mu = 2/trace(R);
        else
            mu = Options.mu;
        end
        if ~isfield(Options,'W0')
            W0 = zeros(M,1);
        else
            W0 = Options.W0;
        end
    end
    % ----->
    W(:,1) = W0; % size is M * 1
    for k = 1:(N-1) % iterate all time stances, in order to caclulate the W of each time
        e = Raw_Signal(k)-Noise_ref(k,:)*W(:,k); % scalar, error in specific time point
        W(:,k+1) = W(:,k) + mu*e*Noise_ref(k,:)'; % calculate the next weights vector (size of next vector is M * 1)
    end
    % now that we have the W, estimate the noise in each instance
    estimatedNoise = zeros(N,1); % preallocate memory
    for i = 1:N
        estimatedNoise(i) = Noise_ref(i,:)*W(:,i);
    end
    Clean_Signal = Raw_Signal - estimatedNoise;  % Subtract the estimated noise (N x 1)
end

