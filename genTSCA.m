function [mapTSCA] = genTSCA(Z)
% Generate TSCA
global fs basis params
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
        [projected,components,~,~,output] = tscaFunc(Z,signal,[noise1 noise2],params.TSCA.gammas,params.TSCA.numProj,params.TSCA.reduceComp); % call TSCA with this signal's parameters
        mapTSCA{i} = postProcess(rshp(components(:,1)));      
    end
    mapTSCA = cat(3,mapTSCA{:});
end

