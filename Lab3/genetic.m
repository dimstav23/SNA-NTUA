function [ list ] = genetic(A,n,p_c,p_m,elitism_n)
%GENETIC 
m=300;
B=cell(300,31);
nei=cell(n);
for i=1:n
    nei{i}=kneighbors(A,i,1);
end
for i=1:300
    for j=1:n
        B{i,1}(j)=nei{j}(randi(size(nei{j},2)));
    end
end
t=1;
prev=0;
prev_t=0;
while t<=30 && prev_t<5 
        D=myfuncell(B,t,n);
        elites=elitism(D);
        for i=1:elitism_n
            B{i,t+1}=B{elites(i),t};
        end
        for i=elitism_n+1:m
            x= rand();
            k=1;
            tmp2=cumsum(D);
            while k<m && x< tmp2(k)/tmp2(m)
                k=k+1;
            end
            B{i,t+1}=B{k,t};
        end
        for i=1:2:299
            if rand()<=p_c
                pos=randi(n-1);
                for k=pos+1:n
                    aux=B{i,t+1}(k);
                    B{i,t+1}(k)=B{i+1,t+1}(k);
                    B{i+1,t+1}(k)=aux;
                end
            end
        end
        for i=1:300
            for k=1:n
                if rand()<p_m
                    B{i,t+1}(k)=nei{k}(randi(size(nei{k},2)));
                end
            end
        end
        max(D)
        if(max(D) ==prev)
            prev_t=prev_t+1;
        else
            prev_t=0;
            prev=max(D);
        end
        t=t+1;
end
[~,j] = max(myfuncell(B,t-1,n));
 Ap=zeros(n);
    for i=1:n
        Ap(B{j,t-1}(i),i)=1;
        Ap(i,B{j,t-1}(i))=1;
    end
    list=find_conn_comp(Ap);
end
function [C] = elitism(B)
    [~,C]=sort(B,'descend');
end
function [C] = myfuncell(B,t,n)
    C=zeros(300,1);
    parfor i=1:300
         C(i)= fittness(B{i,t},n);
    end
end
