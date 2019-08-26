
function  B  = cellcluster( A , n)
%CLUSTERCELL Summary of this function goes here
%   Detailed explanation goes here
B=zeros(n,1)';
b=size(A,2);
for i=1:b
    D=A{i};
    d=size(D,2);
    for j=1:d
        B(D(j))= i;
    end
end
B=B';
end


