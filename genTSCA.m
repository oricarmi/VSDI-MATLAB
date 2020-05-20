function [mapTSCA] = genTSCA(Z)
% Generate TSCA
global fs basis params brn cfn
    noise1.time = eye(params.experiment.T)/params.experiment.T; % autocorrelation matrix of white noise
    if ~isempty(params.TSCA.Noise.freqs) % build oscillatory noises matrices
        for i=1:length(params.TSCA.Noise.freqs)
            noise2(i).time = createToeplitz(params.TSCA.Noise.freqs(i),params.TSCA.Noise.bw,params.TSCA.Noise.numharmonics,params.TSCA.Noise.harmWeights,params.experiment.T);
        end
    else
        noise2 = [];
    end
    mapTSCA = cell(params.experiment.N,1);
    for i=1:params.experiment.N
        signal.time = params.experiment.theoreticalSigs(i,:); % this stimulus theoretical signal
        [projected,components,~,~,output] = tscaFunc(Z,signal,[noise1],params.TSCA.gammas,params.TSCA.numProj,params.TSCA.reduceComp); % call TSCA with this signal's parameters
        [~,I] = max(corr(abs(projected(1:9,:)'),signal.time')); % get the index of component with highest correlation to original signal
        mapTSCA{i} = postProcess(rshp(components(:,I)));      
    end
    mapTSCA = cat(3,mapTSCA{:});
end

