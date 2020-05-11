function [outputArg1,outputArg2] = ScoreAllMaps(mapCorr,mapTSCA,mapT,mapAOF,mapGLM,row,col)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
[~,indCorr] = sort(squeeze(mapCorr(row,col,:))','desc')
[~,indTSCA] = sort(squeeze(mapTSCA(row,col,:))','desc')
[~,indT] = sort(squeeze(mapT(row,col,:))','desc')
[~,indAOF] = sort(squeeze(mapAOF(row,col,:))','desc')
[~,indGLM] = sort(squeeze(mapGLM(row,col,:))','desc')
end

