function [mapAOF] = AvgOfFrms(ZZ)
global params
% Average of frames
mapAOF = cell(params.experiment.N,1);
for i=1:params.experiment.N
    mapAOF{i} = postProcess(mean(ZZ(:,:,(i-1)*params.experiment.T1+params.AOF.numFramesFrom:(i-1)*params.experiment.T1+params.AOF.numFramesUntil),3));
end
mapAOF = cat(3,mapAOF{:});
