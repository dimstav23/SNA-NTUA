clc;clear all;close all;
%%Askisi 1
%problem generated with optimization tool and saved to 3 .mat files
%save('optimproblem.mat','optimproblem');
%save('optimresults.mat','optimresults');
%save('options.mat','options');
%save('optimproblem_double.mat','optimproblem_double');
%save('optimresults_double.mat','optimresults_double');
%save('options_double.mat','options_double');
load('optimproblem.mat');
load('optimresults.mat');
load('options.mat');
load('optimproblem_double.mat');
load('optimresults_double.mat');
load('options_double.mat');
addpath(genpath('../Lab2/untitled folder'));
%addpath(genpath('../Lab1'));
addpath(genpath('./epidemics_code'));
addpath(genpath('/Athanasiou'));
%{
b=0;
c=0;
d=0;
for i=10:10:100
    c=c+1;
    for j=0.3:0.1:0.9
        d=d+1;
        for y=0.01:0.01:0.2
            b=b+1;
            optimproblem.options.PopulationSize = i;
            optimproblem.options.CrossoverFraction = j;
            optimproblem.options.MutationFcn{2} = y;
            t = ga(optimproblem);
            A(b,d,c)=sum(t);

            optimproblem_double.options.PopulationSize = i;
            optimproblem_double.options.CrossoverFraction = j;
            optimproblem_double.options.MutationFcn{2} = y;
            t_double = ga(optimproblem_double);
            A_double(b,d,c)=sum(t_double);
            % print(a);
            %break;
        end
        b=0;
        %break;
    end
    d=0;
end
%}
%%Askisi 2
for i=1:180
    x(i)=rand()*1000; % a random real number between 0-1 * norma . like a pointer
    y(i)=rand()*1000;
end

rgNodes = transpose([x ; y]);

%[football,rggNodeDegree]=rgg(rgNodes,115,250);

%[lesmis,rggNodeDegree]=rgg(rgNodes,77,250);
%[dolphins,rggNodeDegree]=rgg(rgNodes,62,250);
%Lets make the graphs matrix
lesmis =importgml('lesmis.gml');


if isdirected(lesmis)
    B=undir(lesmis,size(lesmis,1));
    lesmis=B;
    clearvars B;
end
%%
maxi=0;
for i=0.7:0.1:0.9
    for j=0.1:0.1:0.2
        for k=1:3
            res=cellcluster(my_ga(lesmis,i,j,k,30,300,5,1),77);
            if QFModul(res,lesmis)>maxi
                maxi=QFModul(res,lesmis);
                maxcluster=res;
            end
        end
    end
end
    
%%            
%res=genetic(lesmis,77,0.7,0.2,2);
figure;
PlotGraph(lesmis,maxcluster);
title('Genetic algorithm for lesmis');
saveas(gcf,'Genetic_lesmis.png');
