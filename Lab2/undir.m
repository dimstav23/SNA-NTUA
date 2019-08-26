function B = undir(A,n)
B=zeros(n);
for i=1:n
    for j=1:n
        if(A(i,j)==1)
            B(i,j)=1;
            B(j,i)=1;
        end
end
end
