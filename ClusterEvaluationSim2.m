function [resDB,resSil] = ClusterEvaluationSim2(maxind)
% Perform cluster similiary analysis from simulated data
U = unique(maxind);
X = cell(1,length(U)); % unique signals detected
for k=2:length(U) % iterate the maps of this method
    [row,col] = find(maxind==U(k));
    X{k-1} = [row,col];
end
[R,db,dbsI,dbns,dbnsI,DBI] = ClusterSimilarity(X(~cellfun(@isempty,X)));
resDB = struct('R',R,'db',db,'dbsI',dbsI,'dbns',dbns,'dbnsI',dbnsI,'DBI',DBI);
figure(1234);
try
    [s,h] = silhouette(cat(1,X{:}),makeClustVector(cellfun(@(x) size(x,1),X))');
    resSil = struct('s',s,'graph',h);
catch
    resSil = struct('s',[],'graph',[]);
end
end

