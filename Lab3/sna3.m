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
%addpath(genpath('/Athanasiou'));

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
            [x,fval,exitflag,output,population,scores] = ga(optimproblem);
            A1(b,d,c)=sum(-1*scores/i);

            %optimproblem_double.options.PopulationSize = i;
            %optimproblem_double.options.CrossoverFraction = j;
            %optimproblem_double.options.MutationFcn{2} = y;
            %t_double = ga(optimproblem_double);
            %A_double(b,d,c)=sum(t_double);
            % print(a);
            %break;
        end
        b=0;
        %break;
    end
    d=0;
end

figure();
plot(linspace(0.01,0.2,20),squeeze(A1(:,3,1))); %crossover 0.5 pupulation 10 
xlabel('Mutation')
ylabel('Scores')
title('crossover = 0.5 , population = 10')
figure();
plot(linspace(0.01,0.2,20),squeeze(A1(:,3,10))); %crossover 0.5 pupulation 100
xlabel('Mutation')
ylabel('Scores')
title('crossover = 0.5 , population = 100')

figure();
plot(linspace(0.3,0.9,7),squeeze(A1(10,:,1))); %mutation 0.1 pupulation 10 
xlabel('Crossover')
ylabel('Scores')
title('mutation = 0.1 , population = 10')
figure();
plot(linspace(0.3,0.9,7),squeeze(A1(10,:,10))); %mutation 0.1 pupulation 100
xlabel('Crossover')
ylabel('Scores')
title('mutation = 0.1 , population = 100')

figure();
plot(linspace(10,100,10),squeeze(A1(1,7,:))); %mutation 0.01 crossover 0.9
xlabel('Population')
ylabel('Scores')
title('mutation = 0.01 , crossover = 0.9')


%% LESMIS
lesmis =importgml('lesmis.gml');
disp('*******Lesmis*********')
if isdirected(lesmis)
    B=undir(lesmis,size(lesmis,1));
    lesmis=B;
    clearvars B;
end
maxi_lesm=0;
for i=0.7:0.1:0.9
    for j=0.1:0.1:0.2
        for k=1:3
            res=cellcluster(my_ga(lesmis,i,j,k,30,300,5,1),77);
            if QFModul(res,lesmis)>maxi_lesm
                maxi_lesm=QFModul(res,lesmis);
                maxcluster_lesm=res;
            end
%            break;
        end
%        break;
    end
%    break;
end           
%res=genetic(lesmis,77,0.7,0.2,2);
PlotGraph(lesmis,maxcluster_lesm);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1); 
X = ['Best Genetic Algorithm Results for Lesmiserables with Modularity Value:' num2str(maxi_lesm)];
title(X);
saveas(gcf,'Genetic_lesmis.png');

%% FOOTBALL
football =importgml('football.gml');
disp('*******Football*********')
if isdirected(football)
    B=undir(football,size(football,1));
    football=B;
    clearvars B;
end
maxi_foot=0;
for i=0.7:0.1:0.9
    for j=0.1:0.1:0.2
        for k=1:3
            res=cellcluster(my_ga(football,i,j,k,30,300,5,1),115);
            if QFModul(res,football)>maxi_foot
                maxi_foot=QFModul(res,football);
                maxcluster_foot=res;
            end
        end
    end
end           
%res=genetic(lesmis,77,0.7,0.2,2);
%figure;
PlotGraph(football,maxcluster_foot);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1); 
X = ['Best Genetic Algorithm Results for Football with Modularity Value:' num2str(maxi_foot)];
title(X);
saveas(gcf,'Genetic_football.png');

%% DOLPHINS
dolphins =importgml('dolphins.gml');
disp('*******Dolphins*********')

if isdirected(dolphins)
    B=undir(dolphins,size(dolphins,1));
    dolphins=B;
    clearvars B;
end
maxi_dol=0;
for i=0.7:0.1:0.9
    for j=0.1:0.1:0.2
        for k=1:3
            res=cellcluster(my_ga(dolphins,i,j,k,30,300,5,1),62);
            if QFModul(res,dolphins)>maxi_dol
                maxi_dol=QFModul(res,dolphins);
                maxcluster_dol=res;
            end
        end
    end
end           
figure;
PlotGraph(dolphins,maxcluster_dol);
frame_h = get(handle(gcf),'JavaFrame');
set(frame_h,'Maximized',1); 
X = ['Best Genetic Algorithm Results for Dolphins with Modularity Value:' num2str(maxi_dol)];
title(X);
saveas(gcf,'Genetic_dolphins.png');

%%SIR
mkdir('SIR');
beta = 10^(-3);
gamma = [10^(-6) 10^(-5) 10^(-4) 10^(-3) 10^(-2) 10^(-1)];
S_0 = 570;
I_0 = 17;
R_0 = 0;
global ind;
for ind=1:6
	SIR(beta,gamma(ind),S_0,I_0,R_0);
end
close all;

%%SIS
clearvars beta
mkdir('SIS');
global alpha N ;
N=1;
for alpha=6:-1:1
	SIS;
end
close all;


