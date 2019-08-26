function Res = MyEgoCent (arr,n)
Res = zeros(n,1);
for i=1:n
    T = [i,kneighbors(arr,i,1)];
    A = subgraph(arr,T);
    ego_A = A*A.*(ones(size(A,1),size(A,2))-A);
    ego_A = triu(ego_A,1);
    ego_A = ego_A .* (1-A);
    ego_A = ego_A(ego_A~=0);
    ego_A = 1./ego_A;
    Res(i) = sum(sum(ego_A));
end
end
