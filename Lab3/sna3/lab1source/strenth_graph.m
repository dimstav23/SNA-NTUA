function [node_str,distri,aver,maxi ] = strenth_graph( A,N)
%STRENTH_GRAPH Summary of this function goes here
%   Detailed explanation goes here
[node_str,~,~]  = degrees(A);
[maxi,distri,~] = cumulativedist(node_str,N);
aver = mean(node_str);
end

