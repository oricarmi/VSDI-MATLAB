function [offsets] = calcOffsets(fitResult)
% calculate the offsets in the brain of the different conditions
% location 1 is (0,0)  
    centers = cellfun(@(x) x(5:6),fitResult(2:end),'UniformOutput',false);
    RefLoc = centers{1};
    offsets = cell2mat(cellfun(@(x) x-RefLoc,centers,'UniformOutput',false));
end

