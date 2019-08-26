function [ Ego_cent ] = ego_cent( matrix,n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for i = 1:n
    V =find(matrix(i,:) ==1);
    S = [i V];
    Adj_temp = subgraph(matrix,S);
    Cent = brandesBetwCentr(Adj_temp);
    Ego_cent(i) = Cent(1);
end

end

