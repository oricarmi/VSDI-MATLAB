function [row2,col2] = convUmp2Num(row,col)
% convert ump to row
global ump 
row2 = row*1000/ump;
col2 = col*1000/ump;
end

