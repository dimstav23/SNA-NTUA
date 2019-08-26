regA=smallw(180,2,0);
[x,y]= getNodeCoordinates(180);
regNodes=[x y];
figure;
gplot(regA,regNodes);
%%

rgA=erdrey(180,750);
figure;
gplot(rgA,regNodes);
%% 
clearvars x y;
%generates random points in LxL
for i=1:180
    x(i)=rand()*1000; % a random real number between 0-1 * norma . like a pointer
    y(i)=rand()*1000;
end

rgNodes = transpose([x ; y]);

[rggA,rggNodeDegree]=rgg(rgNodes,180,250);
figure;
gplot(rggA,rgNodes);
%%

sf=pref(180,4);
figure;
gplot(sf,regNodes);
%%
sw=smallw(180,2,0.3); % the difference with the first is that now we have set the propability.
figure;
gplot(sw,regNodes);
%%

%end of exercise 1.
%Start of exercise 2.
% REG
figure;
s(1) = subplot(2,1,1);
s(2) = subplot(2,1,2);
degreg= degrees(full(regA));
[maximum,cumdistr,degdistr] = cumulativedist(degreg,180);
stem(s(1),[1:maximum],degdistr,'filled','r');
stem(s(2),[1:maximum],cumdistr,'filled','r');
title('Degree and Cumulative distribution for REG');

%   average and variance of degrees
averagedegREG=mean(degreg);
varianceREG=var(degreg);
clearvars s degreg;

%%
%Forfigure;
%RG now
clearvars x y;
s(1) = subplot(2,1,1);
s(2) = subplot(2,1,2);
degreg= degrees(full(rgA));
[maximum,cumdistr,degdistr] = cumulativedist(degreg,180);
stem(s(1),[1:maximum],degdistr,'filled','r');
stem(s(2),[1:maximum],cumdistr,'filled','r');
title('Degree and Cumulative distribution for RG');

%   average and variance of degrees
averagedegRG=mean(degreg);
varianceRG=var(degreg); 
clearvars s degreg;

%% Now for RGG
figure;
s(1) = subplot(2,1,1);
s(2) = subplot(2,1,2);
degreg= degrees(full(rggA));
[maximum,cumdistr,degdistr] = cumulativedist(degreg,180);
stem(s(1),[1:maximum],degdistr,'filled','r');
stem(s(2),[1:maximum],cumdistr,'filled','r');
title('Degree and Cumulative distribution for RGG');

%   average and variance of degrees
averagedegRGG=mean(degreg);
varianceRGG=var(degreg);
clearvars s degreg;
%% SF
figure;
s(1) = subplot(2,1,1);
s(2) = subplot(2,1,2);
degreg= degrees(full(sf));
[maximum,cumdistr,degdistr] = cumulativedist(degreg,180);
stem(s(1),[1:maximum],degdistr,'filled','r');
stem(s(2),[1:maximum],cumdistr,'filled','r');
title('Degree and Cumulative distribution for SF');

%   average and variance of degrees
averagedegSF=mean(degreg);
varianceSF=var(degreg);
clearvars s degreg;

%% SW
figure;
s(1) = subplot(2,1,1);
s(2) = subplot(2,1,2);
degreg= degrees(full(sw));
[maximum,cumdistr,degdistr] = cumulativedist(degreg,180);
stem(s(1),[1:maximum],degdistr,'filled','r');
stem(s(2),[1:maximum],cumdistr,'filled','r');
title('Degree and Cumulative distribution for SW');

%   average and variance of degrees
averagedegSW=mean(degreg);
varianceSW=var(degreg);
clearvars s degreg;

%% EX 3
%reg
W=rand(180)*9 +1 % becouse we need 1 to 10 thats why +1
for i=1:180
    for j=1:180
        W(i,j)=W(j,i);
    end
end
W_reg=zeros(180);
for i=1:180
    for j=1:180
        W_reg(i,j)=W(i,j)*regA(i,j);
    end
end
[reg_node_str,reg_str_dist,reg_avg,reg_max_str]=strenth_graph(full(W_reg),180);
stem([1:reg_max_str],reg_str_dist,'filled','r');
title('Cumulative Strength Distribution-reg');

%% rg

W_rg=zeros(180);
for i=1:180
    for j=1:180
        W_rg(i,j)=W(i,j)*rgA(i,j);
    end
end
[rg_node_str,rg_str_dist,rg_avg,rg_max_str]=strenth_graph(full(W_rg),130);
stem([1:rg_max_str],rg_str_dist,'filled','r');
title('Cumulative Strength Distribution-rg');
%% rgg
W_rgg=zeros(180);
for i=1:180
    for j=1:180
        W_rgg(i,j)=W(i,j)*rggA(i,j);
    end
end
[rgg_node_str,rgg_str_dist,rgg_avg,rgg_max_str]=strenth_graph(full(W_rgg),180);
stem([1:rgg_max_str],rgg_str_dist,'filled','r');
title('Cumulative Strength Distribution-rgg');
%% sf

W_sf=zeros(180);
for i=1:180
    for j=1:180
        W_sf(i,j)=W(i,j)*sf(i,j);
    end
end
[sf_node_str,sf_str_dist,sf_avg,sf_max_str]=strenth_graph(full(W_sf),180);
stem([1:sf_max_str],sf_str_dist,'filled','r');
title('Cumulative Strength Distribution-sf');
%% SW
W_sw=zeros(180);
for i=1:180
    for j=1:180
        W_sw(i,j)=W(i,j)*sw(i,j);
    end
end
[sw_node_str,sw_str_dist,sw_avg,sw_max_str]=strenth_graph(full(W_sw),180);
stem([1:sw_max_str],sw_str_dist,'filled','r');
title('Cumulative Strength Distribution-sw');
%graphs(1)=sym(full(regA));
%graphs(2)=(full(rgA));
%graphs(3)=(full(rggA));
%graphs(4)=(full(sf));
%graphs(5)=(full(sw));
%names(1)='REG';
%names(2)='RG';
%names(3)='RGG';
%names(4)='SF';
%names(5)='SW';
%% EX 4
[path_mean_reg,path_var_reg]=mean_var_path_length(full(regA));
[path_mean_rg,path_var_rg]=mean_var_path_length(full(rgA));
[path_mean_rgg,path_var_rgg]=mean_var_path_length(full(rggA));
[path_mean_sf,path_var_sf]=mean_var_path_length(full(sf));
[path_mean_sw,path_var_sw]=mean_var_path_length(full(sw));
%Kapoai parousiasi na prsothesw...
%% EX5 E2

% REG
figure;
%s(1) = subplot(1,1,1)
[C1,C2,C]= clust_coeff(full(regA));
%stem(s(1),[1:180],C,'filled','r');
local_clust_avg_reg=C2;
cdfplot(C);
title('Clustering Coefficient distribution for REG');
hold on;
clearvars C1 C2 C;
%% RG
figure;
%s(1) = subplot(1,1,1);
[C1,C2,C]= clust_coeff(full(rgA));
%stem(s(1),[1:180],C,'filled','r');
cdfplot(C);
title('Clustering Coefficient distribution for RG');
hold on;
local_clust_avg_rg=C2;
clearvars C1 C2 C;
%% RGG
figure;
[C1,C2,C]= clust_coeff(full(rggA));
%stem(s(1),[1:180],C,'filled','r');
cdfplot(C);
title('Clustering Coefficient distribution for RGG');
hold on;
local_clust_avg_rgg=C2;
clearvars C1 C2 C;
%% SF
figure;
[C1,C2,C]= clust_coeff(full(sf));
cdfplot(C);
title('Clustering Coefficient distribution for SF');
hold on;
local_clust_avg_sf=C2;
clearvars C1 C2 C;
%% SW
figure;
[C1,C2,C]= clust_coeff(full(sw));
cdfplot(C);
title('Clustering Coefficient distribution for SW');
hold on;
local_clust_avg_sw=C2;
clearvars C1 C2 C;

%% ex Z
% REG
figure;

[degree_cent_reg,~,~]=degrees(regA);
closs_cent_reg=closeness(regA);
betw_cent_reg=node_betweenness_faster(regA);
eigen_cent_reg=eigencentrality(full(regA));
[~,degree_cent_reg1]=cumulativecentrality(degree_cent_reg,130);
[~,betw_cent_reg1]=cumulativecentrality(betw_cent_reg,130);
[~,closs_cent_reg1]=cumulativecentrality(closs_cent_reg,130);
[~,eigen_cent_reg1]=cumulativecentrality(eigen_cent_reg,130);
subplot(4,1,1);
cdfplot(full(degree_cent_reg1));
title('Deegree centrality REG');
subplot(4,1,2);
cdfplot(betw_cent_reg1);
title('Closeness centrality REG');
subplot(4,1,3);
cdfplot(closs_cent_reg1);
title('Betwwness centrality REG');
subplot(4,1,4);
cdfplot(full(eigen_cent_reg1));
title('Eigenvector centality REG');
hold all;
degree_mean_cent_reg=mean(degree_cent_reg);
closs_avg_cent_reg=mean(closs_cent_reg);
bwtw_avg_cent_reg=mean(betw_cent_reg);
eigen_avg_cent_reg=mean(eigen_cent_reg);
%% RG
figure;
%s(1) = subplot(4,1,1);
%s(2) = subplot(4,1,2);
%s(3) = subplot(4,1,3);
%s(4) = subplot(4,1,4);
%[degree_cent_rg,~,~]=degrees(rgA);
%closs_cent_rg=closeness(rgA);
%betw_cent_rg=node_betweenness_faster(rgA);
%eigen_cent_rg=eigencentrality(full(rgA));
%[~,degree_cent_rg1]=cumulativecentrality(degree_cent_rg,180);
%[~,betw_cent_rg1]=cumulativecentrality(betw_cent_rg,180);
%[~,closs_cent_rg1]=cumulativecentrality(closs_cent_rg,180);
%[~,eigen_cent_rg1]=cumulativecentrality(eigen_cent_rg,180);
%stem(s(1),1:101,full(degree_cent_rg1).','filled','r');
%stem(s(2),1:101,closs_cent_rg1,'filled','r');
%stem(s(3),1:101,betw_cent_rg1,'filled','r');
%stem(s(4),1:101,full(eigen_cent_rg1).','filled','r');
subplot(4,1,1);
cdfplot(full(degree_cent_rg));
title('Deegree centrality RG');
subplot(4,1,2);
cdfplot(betw_cent_rg);
title('Closeness centrality RG');
subplot(4,1,3);
cdfplot(closs_cent_rg);
title('Betwwness centrality RG');
subplot(4,1,4);
cdfplot(full(eigen_cent_rg));
title('Eigenvector centality RG');
hold all;
degree_mean_cent_rg=mean(degree_cent_rg);
closs_avg_cent_rg=mean(closs_cent_rg);
bwtw_avg_cent_rg=mean(betw_cent_rg);
eigen_avg_cent_rg=mean(eigen_cent_rg);

%% RGG
figure;

[degree_cent_rgg,~,~]=degrees(rggA);
closs_cent_rgg=closeness(rggA);
betw_cent_rgg=node_betweenness_faster(rggA);
eigen_cent_rgg=eigencentrality(full(rggA));
[~,degree_cent_rgg1]=cumulativecentrality(degree_cent_rgg,180);
[~,betw_cent_rgg1]=cumulativecentrality(betw_cent_rgg,180);
[~,closs_cent_rgg1]=cumulativecentrality(closs_cent_rgg,180);
[~,eigen_cent_rgg1]=cumulativecentrality(eigen_cent_rgg,180);
subplot(4,1,1);
cdfplot(full(degree_cent_rgg));
title('Deegree centrality RGG');
subplot(4,1,2);
cdfplot(betw_cent_rgg);
title('Closeness centrality RGG');
subplot(4,1,3);
cdfplot(closs_cent_rgg);
title('Betwwness centrality RGG');
subplot(4,1,4);
cdfplot(full(eigen_cent_rgg));
title('Eigenvector centality RGG');
hold all;
degree_mean_cent_rgg=mean(degree_cent_rgg);
closs_avg_cent_rgg=mean(closs_cent_rgg);
bwtw_avg_cent_rgg=mean(betw_cent_rgg);
eigen_avg_cent_rgg=mean(eigen_cent_rgg);
%% SF
figure;
s(1) = subplot(4,1,1);
s(2) = subplot(4,1,2);
s(3) = subplot(4,1,3);
s(4) = subplot(4,1,4);
[degree_cent_sf,~,~]=degrees(sf);
closs_cent_sf=closeness(sf);
betw_cent_sf=node_betweenness_faster(sf);
eigen_cent_sf=eigencentrality(full(sf));
[~,degree_cent_sf1]=cumulativecentrality(degree_cent_sf,180);
[~,betw_cent_sf1]=cumulativecentrality(betw_cent_sf,180);
[~,closs_cent_sf1]=cumulativecentrality(closs_cent_sf,180);
[~,eigen_cent_sf1]=cumulativecentrality(eigen_cent_sf,180);
subplot(4,1,1);
cdfplot(full(degree_cent_sf));
title('Deegree centrality SF');
subplot(4,1,2);
cdfplot(betw_cent_sf);
title('Closeness centrality SF');
subplot(4,1,3);
cdfplot(closs_cent_sf);
title('Betwwness centrality SF');
subplot(4,1,4);
cdfplot(full(eigen_cent_sf));
title('Eigenvector centality SF');
hold all;
degree_mean_cent_sf=mean(degree_cent_sf);
closs_avg_cent_sf=mean(closs_cent_sf);
bwtw_avg_cent_sf=mean(betw_cent_sf);
eigen_avg_cent_sf=mean(eigen_cent_sf);
%% SW
figure;
s(1) = subplot(4,1,1);
s(2) = subplot(4,1,2);
s(3) = subplot(4,1,3);
s(4) = subplot(4,1,4);
[degree_cent_sw,~,~]=degrees(sw);
closs_cent_sw=closeness(sw);
betw_cent_sw=node_betweenness_faster(sw);
eigen_cent_sw=eigencentrality(full(sw));
[~,degree_cent_sw1]=cumulativecentrality(degree_cent_sw,180);
[~,betw_cent_sw1]=cumulativecentrality(betw_cent_sw,180);
[~,closs_cent_sw1]=cumulativecentrality(closs_cent_sw,180);
[~,eigen_cent_sw1]=cumulativecentrality(eigen_cent_sw,180);
subplot(4,1,1);
cdfplot(full(degree_cent_sw));
title('Deegree centrality SW');
subplot(4,1,2);
cdfplot(betw_cent_sw);
title('Closeness centrality SW');
subplot(4,1,3);
cdfplot(closs_cent_sw);
title('Betwwness centrality SW');
subplot(4,1,4);
cdfplot(full(eigen_cent_sw));
title('Eigenvector centality SW');
hold all;
degree_mean_cent_sw=mean(degree_cent_sw);
closs_avg_cent_sw=mean(closs_cent_sw);
bwtw_avg_cent_sw=mean(betw_cent_sw);
eigen_avg_cent_sw=mean(eigen_cent_sw);

%%
j=1;
for i=0:0.1:1
    sw_I=smallw(180,2,i);
    ave_path_sw_I(j)=ave_path_length(sw_I);
    [~,average_clust_sw_I(j),~]=clust_coeff(sw_I);
    j=j+1;
end
%%
regA11=smallw(180,2,1);
[x,y]= getNodeCoordinates(180);
regNodes11=[x y];
figure;
gplot(regA11,regNodes11);
