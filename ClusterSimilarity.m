function [R,dbs,dbsI,dbns,dbnsI,DBI] = ClusterSimilarity(X,qp)
% Calculate the Davies-Bouldin similarity index
% input: X - cell array of clusters. each cell (cluster) - row=observations, cols= dimensions
% for 2D: 1st column - rows indices of cluster, 2nd column - col indices of cluster
    switch nargin
        case 1
            q=2;p=2; % default p,q values
        case 2
            if isempty(qp(1))
                q = 2;
            end
            if isempty(qp(2))
                p=2;
            end
        otherwise
            disp('not enough input arguments');
            return
    end
    k = length(X);
    A = zeros(k,2); % centroids (rows - clusters, cols - dimensions)
    S = zeros(k,1); % scatters (rows - clusters, cols - dimensions)
    M = zeros(k); % distance measures
    for i=1:k
        A(i,:) = calcCentroid(X{i});
        S(i) = calcSi(X{i},A(i,:),size(X{i},1),q);
    end
    for i = 1:k
        for j = i+1:k
            M(i,j) = calcM(A(i,:),A(j,:),p);
            M(j,i) = M(i,j);
        end
    end
    % Davies-Bouldin
    R = zeros(k);
    dbs=zeros(1,k);dbsI = dbs; dbns = dbs; dbnsI = dbs;
    for i = 1:k
        for j = i+1:k
            R(i,j) = (S(i) + S(j))/M(i,j);
            R(j,i) = R(i,j);
        end
        [dbs(i),dbsI(i)] = max(R(i,:));
        [tmp,tmpI] = mink(R(i,:),2);
        dbns(i) = tmp(2); dbnsI(i) = tmpI(2);
    end
    DBI = mean(dbs);% mean Davies-Bouldin
    R = R + eye(size(R,1)); % add ones to diagonal (similarity between a cluster and itself is 1)

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBMETHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Si = calcSi(X,A,T,q)
    Si = (1/T)*(sum(sum(abs(X-A).^q))).^(1/q);
end
function Ai = calcCentroid(X)
    Ai = mean(X);
end
function Mij = calcM(Ai,Aj,p)
    Mij = sum(abs(Ai-Aj).^p).^(1/p);
end