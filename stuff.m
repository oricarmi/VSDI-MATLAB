addpath('C:\Users\orica\Dropbox\fcns_and_decript');
addpath('C:\Users\orica\Dropbox\master degree\codes');
path = 'C:\Users\Ori\Desktop\Ori\2nd degree\matlab codez\vsdi - matlab\comparison results';
path = 'C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\comparison results';
files = dir(path);
global brn lgn ump
brn = zeros(270,327);
for i=7:length(files) % iterate files
    if contains(files(i).name,'181218')
        continue
    end
    load(fullfile(files(i).folder,files(i).name)); % load summary
    result = Summary.result;
    paramz = Summary.params;
    optimalMaps = Summary.params.experiment.optimalMaps.orig;
    N = Summary.params.experiment.N;
    X = cell(N,7); % cell array of clusters (rows - maps, cols - methods)
    fn = fieldnames(result);
    RR = mean(Summary.params.experiment.optimalMaps.orig,3)>prctile(reshape(mean(Summary.params.experiment.optimalMaps.orig,3),[],1),80);
    for j=1:length(fn)+1 % iterate the methods and the optimal maps and perform cluster analysis between the maps
        if j~=length(fn)+1 % if not in optimal maps
            for k=1:N % iterate the maps of this method
                thisMap = result.(fn{j}).maps(:,:,k).*double(RR);
                [row,col] = find(thisMap>prctile(reshape(thisMap,[],1),99));
%                 mapp(:,:,k) = thisMap>prctile(reshape(thisMap,[],1),99);
                X{k,j} = [row,col];
            end
            [R,dbs,dbsI,dbns,dbnsI,DBI] = ClusterSimilarity(X(:,j));
            result.(fn{j}).clusterEval.R = R;
            result.(fn{j}).clusterEval.dbs = dbs;
            result.(fn{j}).clusterEval.dbsI = dbsI;
            result.(fn{j}).clusterEval.dbns = dbns;
            result.(fn{j}).clusterEval.dbnsI = dbnsI;
            result.(fn{j}).clusterEval.DBI = DBI;
        else % it is optimal maps
            for k=1:N % iterate the optimal maps
                [row,col] = find(optimalMaps(:,:,k)>prctile(reshape(optimalMaps(:,:,k),[],1),99));
%                 mapp(:,:,k) = optimalMaps(:,:,k)>prctile(reshape(optimalMaps(:,:,k),[],1),99);
                X{k,j} = [row,col];
            end
            [R,dbs,dbsI,dbns,dbnsI,DBI] = ClusterSimilarity(X(:,j));
            paramz.experiment.optimalMaps.clusterEval.R = R;
            paramz.experiment.optimalMaps.clusterEval.dbs = dbs;
            paramz.experiment.optimalMaps.clusterEval.dbsI = dbsI;
            paramz.experiment.optimalMaps.clusterEval.dbns = dbns;
            paramz.experiment.optimalMaps.clusterEval.dbnsI = dbnsI;
            paramz.experiment.optimalMaps.clusterEval.DBI = DBI;
        end
    end
    for k=1:N % iterate the maps and perform cluster analysis between the methods
        [R,dbs,dbsI,dbns,dbnsI,DBI] = ClusterSimilarity(X(k,:));
        result.clusterEvalAll{k}.R = R;
        result.clusterEvalAll{k}.dbs = dbs;
        result.clusterEvalAll{k}.dbsI = dbsI;
        result.clusterEvalAll{k}.dbns = dbns;
        result.clusterEvalAll{k}.dbnsI = dbnsI;
        result.clusterEvalAll{k}.DBI = DBI;
    end 
    Summary.result = result;
    Summary.params = paramz;
    uisave('Summary',Summary.description)
end

%% calc real snr
ZSig = zeros(m*m,T); % preallocate memory
for i = 1:T
    ZSig(:,i) = reshape( ...
        signals(1).time(i)*signals(1).space+signals(2).time(i)*signals(2).space+...
        signals(3).time(i)*signals(3).space+signals(4).time(i)*signals(4).space+...
        signals(5).time(i)*signals(5).space,[],1);
end
ZNoise = zeros(m*m,T); % preallocate memory
for i = 1:T
    ZNoise(:,i) = reshape( ...
        noises(1).time(i)*noises(1).space+...
        noises(2).time(i)*noises(2).space+...
        noises(3).time(i)*noises(3).space,[],1);
end
%%
addpath('C:\Users\orica\Dropbox\fcns_and_decript');
addpath('C:\Users\orica\Dropbox\master degree\codes');
path = 'C:\Users\Ori\Desktop\Ori\2nd degree\matlab codez\vsdi - matlab\comparison results';
% path = 'C:\Users\orica\OneDrive\Desktop\2nd degree\matlab codez\matlab - vsdi\comparison results';
files = dir(path);
global brn lgn ump
brn = zeros(270,327);
allPerformance = []; allDBI = [];
for i=3:length(files) % iterate files
    load(fullfile(files(i).folder,files(i).name)); % load summary
    result = Summary.result;
