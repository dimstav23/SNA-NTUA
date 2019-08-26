function  B  = clustercell( A , n,p)
%CLUSTERCELL Summary of this function goes here
%   Detailed explanation goes here
B=zeros(n,1)';
b=size(A{p});
for i=1:b(2)
    D=A{p}{i};
    d=size(D)
    for j=1:d(2)
        B(D(j))= i;
    end
end
B=B'
