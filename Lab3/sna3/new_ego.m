function [ Ego_cent ] = new_ego( A,L )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
for j = 1: L
  v =find(A(j,:) ==1);
  s = [j v];
  adj_temp = subgraph(A,s);
  Atemp = (adj_temp^2).*(ones(size(adj_temp)) - adj_temp);
  U = triu(Atemp)-diag(diag(Atemp));    
  Ego_cent(j) = sum(1./U(U~=0));
end


end
