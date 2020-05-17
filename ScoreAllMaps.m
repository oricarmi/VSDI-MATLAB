function [outputMax,outputInd] = ScoreAllMaps(mapTSCA,mapT,mapAOF,mapCorr,mapGLM,mapNadav,row,col)
% get the value of each pixel sorted per method
    if mod(row,1) || mod(col,1) || (col<=7 && row<=7)
        [row,col] = convUmp2Num(row,col);
    end
    [maxTSCA,indTSCA] = sort(squeeze(mapTSCA(row,col,:))','desc');
    [maxTmax,indTmax] = sort(squeeze(mapT(row,col,:))','desc');
    [maxAOF,indAOF] = sort(squeeze(mapAOF(row,col,:))','desc');
    [maxCorr,indCorr] = sort(squeeze(mapCorr(row,col,:))','desc');
    [maxGLM,indGLM] = sort(squeeze(mapGLM(row,col,:))','desc');
    [maxNadav,indNadav] = sort(squeeze(mapNadav(row,col,:))','desc');
    outputInd = [indTSCA;indAOF;indCorr;indGLM;indTmax;indNadav]
    outputMax = [maxTSCA;maxAOF;maxCorr;maxGLM;maxTmax;maxNadav]
end

