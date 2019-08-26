function [node_str,distri,aver,maxi ] = strength_graph( A,N)
%STRENGTH_GRAPH Summary of this function goes here
%   Detailed explanation goes here
[node_str,~,~]  = degrees(A);
[maxi,distri,~] = cumulativedist(node_str,N);
aver = mean(node_str);
end

