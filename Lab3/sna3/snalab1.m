%% Vavouliotis Giorgos (AM:03112083)

clear all; close all; clc;

%% Step A
mkdir StepA

% REG
regA  = smallw(130,2,0); % vathmos d = 4 alla h smallw kanei d = 2*d
[cord1,cord2] = getNodeCoordinates(130);
regNodes = [cord1 cord2];
cd StepA
figure(1);
gplot(regA,regNodes);
title('REG');
print -djpeg REG.jpeg; 
cd ../

% RG(ER)
rgA = erdrey(130,750);
cd StepA
figure(2);
gplot(rgA,regNodes);
title('RG(ER)');
print -djpeg RG.jpeg; 
cd ..

% RGG
clearvars x y;
for i=1:130
    x(i)=rand()*1000; % a random real number between 0-1 * 1000
    y(i)=rand()*1000;
end

rgNodes = [x ; y]' ;
[rggA,rggNodeDegree] = rgg(rgNodes,130,250);
cd StepA
figure(3);
gplot(rggA,rgNodes);
title('RGG');
print -djpeg RGG.jpeg; 
cd ..

%SF(BA)
sf = pref(130,4);
cd StepA
figure(4);
gplot(sf,regNodes);
title('SF(BA)');
print -djpeg SF.jpeg; 
cd ..

%SW(WS)
sw=smallw(130,2,0.3); 
cd StepA
figure(5);
gplot(sw,regNodes);
title('SW(WS)');
print -djpeg SW.jpeg; 
cd ..

close all ; 

%% Step B
mkdir StepB
cd StepB; 
mkdir REG; mkdir RG; mkdir RGG; mkdir SF; mkdir SW;
cd ..

% REG
degreg = degrees(regA);
[supremum,cumdist,degdist] = cumulativedist(degreg,130);
cd StepB/REG;
figure(6);
plot([1:supremum],degdist);
title('Degree distribution for REG');
print -djpeg Deegre.jpeg;
figure(7);
stem([1:supremum],degdist,'filled','r');
title('Degree distribution for REG(discrete)');
print -djpeg DeegreDisc.jpeg;

figure(8);
plot([1:supremum],cumdist);
title('Cumulative distribution for REG');
print -djpeg Cumulative.jpeg;
figure(9);
stem([1:supremum],cumdist,'filled','r');
title('Cumulative distribution for REG(discrete)');
print -djpeg CumulativeDisc.jpeg;

cd ../../

temp = degrees(full(regA));
averagedegREG = mean(temp); % average
varianceREG   = var(temp);  % variance

close all;

% RG 
degrg = degrees(rgA);
[supremum,cumdist,degdist] = cumulativedist(degrg,130);
cd StepB/RG;
figure(10);
plot([1:supremum],degdist);
title('Degree distribution for RG');
print -djpeg Deegre.jpeg;
figure(11);
stem([1:supremum],degdist,'filled','r');
title('Degree distribution for RG(discrete)');
print -djpeg DeegreDisc.jpeg;

figure(12);
plot([1:supremum],cumdist);
title('Cumulative distribution for RG');
print -djpeg Cumulative.jpeg;
figure(13);
stem([1:supremum],cumdist,'filled','r');
title('Cumulative distribution for RG(discrete)');
print -djpeg CumulativeDisc.jpeg;

cd ../../

clearvars temp ;
temp = degrees(full(rgA));
averagedegRG = mean(temp); % average
varianceRG   = var(temp);  % variance


%RGG
degrgg = degrees(rggA);
[supremum,cumdist,degdist] = cumulativedist(degrgg,130);
cd StepB/RGG;
figure(14);
plot([1:supremum],degdist);
title('Degree distribution for RGG');
print -djpeg Deegre.jpeg;
figure(15);
stem([1:supremum],degdist,'filled','r');
title('Degree distribution for RGG(discrete)');
print -djpeg DeegreDisc.jpeg;

figure(16);
plot([1:supremum],cumdist);
title('Cumulative distribution for RGG');
print -djpeg Cumulative.jpeg;
figure(17);
stem([1:supremum],cumdist,'filled','r');
title('Cumulative distribution for RGG(discrete)');
print -djpeg CumulativeDisc.jpeg;

cd ../../

clearvars temp ;
temp = degrees(full(rggA));
averagedegRGG = mean(temp); % average
varianceRGG   = var(temp);  % variance

%SF 
degsf = degrees(sf);
[supremum,cumdist,degdist] = cumulativedist(degsf,130);
cd StepB/SF;
figure(18);
plot([1:supremum],degdist);
title('Degree distribution for SF');
print -djpeg Deegre.jpeg;
figure(19);
stem([1:supremum],degdist,'filled','r');
title('Degree distribution for SF(discrete)');
print -djpeg DeegreDisc.jpeg;

figure(20);
plot([1:supremum],cumdist);
title('Cumulative distribution for SF');
print -djpeg Cumulative.jpeg;
figure(21);
stem([1:supremum],cumdist,'filled','r');
title('Cumulative distribution for SF(discrete)');
print -djpeg CumulativeDisc.jpeg;

cd ../../

clearvars temp ;
temp = degrees(full(sf));
averagedegSF = mean(temp); % average
varianceSF   = var(temp);  % variance

%SW
degsw = degrees(sw);
[supremum,cumdist,degdist] = cumulativedist(degsw,130);
cd StepB/SW;
figure(22);
plot([1:supremum],degdist);
title('Degree distribution for SW');
print -djpeg Deegre.jpeg;
figure(23);
stem([1:supremum],degdist,'filled','r');
title('Degree distribution for SW(discrete)');
print -djpeg DeegreDisc.jpeg;

figure(24);
plot([1:supremum],cumdist);
title('Cumulative distribution for SW');
print -djpeg Cumulative.jpeg;
figure(25);
stem([1:supremum],cumdist,'filled','r');
title('Cumulative distribution for SW(discrete)');
print -djpeg CumulativeDisc.jpeg;

cd ../../

clearvars temp ;
temp = degrees(full(sw));
averagedegSW = mean(temp); % average
varianceSW   = var(temp);  % variance

close all;

%% Step C
mkdir StepC
cd StepC; 
mkdir REG; mkdir RG; mkdir RGG; mkdir SF; mkdir SW;
cd ..

%rand(130) returns an array 130X130 and +1 is because i want weights 1-10 
weightsarray = rand(130)*9 +1 ;
for i=1:130
    for j=1:130
        weightsarray(i,j)=weightsarray(j,i); % make it symmetric
    end
end

% REG
wreg = zeros(130);
for i=1:130
    for j=1:130
        wreg(i,j) = weightsarray(i,j)*regA(i,j);
    end
end

[reg_node_str,reg_str_dist,reg_avg,reg_max_str] = strenth_graph(wreg,130);
cd StepC/REG;
figure(26);
plot([1:reg_max_str],reg_str_dist);
title('Cumulative Strength Distribution-REG');
print -djpeg Strength.jpeg;
figure(27);
stem([1:reg_max_str],reg_str_dist,'filled','r');
title('Cumulative Strength Distribution-REG(discrete)');
print -djpeg StrengthDisc.jpeg;
cd ../..

%RG
wrg = zeros(130);
for i=1:130
    for j=1:130
        wrg(i,j) = weightsarray(i,j)*rgA(i,j);
    end
end

[reg_node_str,reg_str_dist,reg_avg,reg_max_str] = strenth_graph(wrg,130);
cd StepC/RG;
figure(28);
plot([1:reg_max_str],reg_str_dist);
title('Cumulative Strength Distribution-RG');
print -djpeg Strength.jpeg;
figure(29);
stem([1:reg_max_str],reg_str_dist,'filled','r');
title('Cumulative Strength Distribution-RG(discrete)');
print -djpeg StrengthDisc.jpeg;
cd ../..


%RGG
wrgg = zeros(130);
for i=1:130
    for j=1:130
        wrgg(i,j) = weightsarray(i,j)*rggA(i,j);
    end
end

[reg_node_str,reg_str_dist,reg_avg,reg_max_str] = strenth_graph(wrgg,130);
cd StepC/RGG;
figure(30);
plot([1:reg_max_str],reg_str_dist);
title('Cumulative Strength Distribution-RGG');
print -djpeg Strength.jpeg;
figure(31);
stem([1:reg_max_str],reg_str_dist,'filled','r');
title('Cumulative Strength Distribution-RGG(discrete)');
print -djpeg StrengthDisc.jpeg;
cd ../..

%SF
wsf = zeros(130);
for i=1:130
    for j=1:130
        wsf(i,j) = weightsarray(i,j)*sf(i,j);
    end
end

[reg_node_str,reg_str_dist,reg_avg,reg_max_str] = strenth_graph(wsf,130);
cd StepC/SF;
figure(32);
plot([1:reg_max_str],reg_str_dist);
title('Cumulative Strength Distribution-SF');
print -djpeg Strength.jpeg;
figure(33);
stem([1:reg_max_str],reg_str_dist,'filled','r');
title('Cumulative Strength Distribution-SF(discrete)');
print -djpeg StrengthDisc.jpeg;
cd ../..


%SW
wsw = zeros(130);
for i=1:130
    for j=1:130
        wsw(i,j) = weightsarray(i,j)*sw(i,j);
    end
end

[reg_node_str,reg_str_dist,reg_avg,reg_max_str] = strenth_graph(wsw,130);
cd StepC/SW;
figure(34);
plot([1:reg_max_str],reg_str_dist);
title('Cumulative Strength Distribution-SW');
print -djpeg Strength.jpeg;
figure(35);
stem([1:reg_max_str],reg_str_dist,'filled','r');
title('Cumulative Strength Distribution-SW(discrete)');
print -djpeg StrengthDisc.jpeg;
cd ../..

%% Step D
[path_mean_reg,path_var_reg] = ave_path_length(regA);
[path_mean_rg,path_var_rg]   = ave_path_length(rgA);
[path_mean_rgg,path_var_rgg] = ave_path_length(rggA);
[path_mean_sf,path_var_sf]   = ave_path_length(sf);
[path_mean_sw,path_var_sw]   = ave_path_length(sw);

AllInOne      = zeros(5,2);
AllInOne(1,1) = path_mean_reg; AllInOne(1,2) = path_var_reg; 
AllInOne(2,1) = path_mean_rg;  AllInOne(2,2) = path_var_rg; 
AllInOne(3,1) = path_mean_rgg; AllInOne(3,2) = path_var_rgg; 
AllInOne(4,1) = path_mean_sf; AllInOne(4,2) = path_var_sf;
AllInOne(5,1) = path_mean_sw; AllInOne(5,2) = path_var_sw;

%% Step E
mkdir StepE
cd StepE; 
mkdir REG; mkdir RG; mkdir RGG; mkdir SF; mkdir SW;
cd ..

% REG
[C1,C2,C]= clust_coeff(full(regA));
local_clust_avg_reg=C2;
figure(36);
cd StepE/REG
cdfplot(C);
title('Clustering Coefficient distribution for REG');
print -djpeg CCREG.jpeg;
clearvars C1 C2 C;
cd ../..

% RG
[C1,C2,C]= clust_coeff(full(rgA));
local_clust_avg_rg=C2;
cd StepE/RG
figure(37);
cdfplot(C);
title('Clustering Coefficient distribution for RG');
print -djpeg CCRG.jpeg;
clearvars C1 C2 C;
cd ../..

% RGG
[C1,C2,C]= clust_coeff(full(rggA));
local_clust_avg_rgg=C2;
cd StepE/RGG
figure(38);
cdfplot(C);
title('Clustering Coefficient distribution for RGG');
print -djpeg CCRGG.jpeg;
clearvars C1 C2 C;
cd ../..

% SF
[C1,C2,C]= clust_coeff(full(sf));
local_clust_avg_sf=C2;
cd StepE/SF
figure(39);
cdfplot(C);
title('Clustering Coefficient distribution for SF');
print -djpeg CCSF.jpeg;
clearvars C1 C2 C;
cd ../..

% SW
[C1,C2,C]= clust_coeff(full(sw));
local_clust_avg_sw=C2;
cd StepE/SW
figure(40);
cdfplot(C);
title('Clustering Coefficient distribution for SW');
print -djpeg CCSW.jpeg;
clearvars C1 C2 C;
cd ../..

%% Step Z
mkdir StepZ
cd StepZ; 
mkdir REG; mkdir RG; mkdir RGG; mkdir SF; mkdir SW;
cd ..

% REG
[degree_cent_reg,~,~] = degrees(regA);
[~,degree_cent_reg1]  = cumulativecentrality(degree_cent_reg,130);

closs_cent_reg        = closeness(regA);
[~,closs_cent_reg1]   = cumulativecentrality(closs_cent_reg,130);

betw_cent_reg         = node_betweenness_faster(regA);
[~,betw_cent_reg1]    = cumulativecentrality(betw_cent_reg,130);

eigen_cent_reg = eigencentrality(full(regA));
[~,eigen_cent_reg1]  = cumulativecentrality(eigen_cent_reg,130);

figure(41);
hold on;
subplot(2,1,1); cdfplot(degree_cent_reg1); title('Deegree centrality REG');
subplot(2,1,2); cdfplot(betw_cent_reg1);   title('Closeness centrality REG');
hold off;
cd StepZ/REG 
print -djpeg CentralityReg.jpeg;
cd ../..

figure(42);
hold on;
subplot(2,1,1); stem([1:101],degree_cent_reg1,'filled','b'); title('Deegree centrality REG');
subplot(2,1,2); stem([1:101],betw_cent_reg1,'filled','b');  title('Closeness centrality REG');
hold off;
cd StepZ/REG 
print -djpeg CentralityRegStem.jpeg;
cd ../..

figure(43);
hold on;
subplot(2,1,1); cdfplot(closs_cent_reg1);  title('Betwwness centrality REG');
subplot(2,1,2); cdfplot(eigen_cent_reg1);  title('Eigenvector centality REG');
hold off;
cd StepZ/REG 
print -djpeg CentralityReg1.jpeg;
cd ../..

figure(44);
hold on;
subplot(2,1,1); stem([1:101],closs_cent_reg1,'filled','b');  title('Betwwness centrality REG');
subplot(2,1,2); stem([1:101],eigen_cent_reg1,'filled','b');  title('Eigenvector centality REG');
hold off;
cd StepZ/REG 
print -djpeg CentralityReg1Stem.jpeg;
cd ../..

degree_mean_cent_reg=mean(degree_cent_reg);
closs_avg_cent_reg=mean(closs_cent_reg);
bwtw_avg_cent_reg=mean(betw_cent_reg);
eigen_avg_cent_reg=mean(eigen_cent_reg);


% RG
[degree_cent_rg,~,~] = degrees(full(rgA));
[~,degree_cent_rg1]  = cumulativecentrality(degree_cent_rg,130);
stem([1:101],degree_cent_rg1,'filled','b');
closs_cent_rg = closeness(rgA);
[~,closs_cent_rg1]  = cumulativecentrality(closs_cent_rg,130);

betw_cent_rg  = node_betweenness_faster(rgA);
[~,betw_cent_rg1]   = cumulativecentrality(betw_cent_rg,130);

eigen_cent_rg = eigencentrality(full(rgA));
[~,eigen_cent_rg1]  = cumulativecentrality(eigen_cent_rg,130);

figure(45);
hold on;
subplot(2,1,1); cdfplot(degree_cent_rg1); title('Deegree centrality RG');
subplot(2,1,2); cdfplot(betw_cent_rg1);   title('Closeness centrality RG');
hold off;
cd StepZ/RG 
print -djpeg CentralityRg.jpeg;
cd ../..

figure(46);
hold on;
subplot(2,1,1); cdfplot(closs_cent_rg1);  title('Betwwness centrality RG');
subplot(2,1,2); cdfplot(eigen_cent_rg1);  title('Eigenvector centality RG');
hold off;
cd StepZ/RG 
print -djpeg CentralityRg1.jpeg;
cd ../..

degree_mean_cent_rg=mean(degree_cent_rg);
closs_avg_cent_rg=mean(closs_cent_rg);
bwtw_avg_cent_rg=mean(betw_cent_rg);
eigen_avg_cent_rg=mean(eigen_cent_rg);


% RGG
[degree_cent_rgg,~,~] = degrees(full(rggA));
[~,degree_cent_rgg1]  = cumulativecentrality(degree_cent_rgg,130);

closs_cent_rgg = closeness(rggA);
[~,closs_cent_rgg1]  = cumulativecentrality(closs_cent_rgg,130);

betw_cent_rgg  = node_betweenness_faster(rggA);
[~,betw_cent_rgg1]   = cumulativecentrality(betw_cent_rgg,130);

eigen_cent_rgg = eigencentrality(full(rggA));
[~,eigen_cent_rgg1]  = cumulativecentrality(eigen_cent_rgg,130);

figure(47);
hold on;
subplot(2,1,1); cdfplot(degree_cent_rg1); title('Deegree centrality RGG');
subplot(2,1,2); cdfplot(betw_cent_rg1);   title('Closeness centrality RGG');
hold off;
cd StepZ/RGG 
print -djpeg CentralityRgg.jpeg;
cd ../..

figure(48);
hold on;
subplot(2,1,1); cdfplot(closs_cent_rg1);  title('Betwwness centrality RGG');
subplot(2,1,2); cdfplot(eigen_cent_rg1);  title('Eigenvector centality RGG');
hold off;
cd StepZ/RGG 
print -djpeg CentralityRgg1.jpeg;
cd ../..

degree_mean_cent_rgg=mean(degree_cent_rgg);
closs_avg_cent_rgg=mean(closs_cent_rgg);
bwtw_avg_cent_rgg=mean(betw_cent_rgg);
eigen_avg_cent_rgg=mean(eigen_cent_rgg);


% SF
[degree_cent_sf,~,~] = degrees(full(sf));
[~,degree_cent_sf1]  = cumulativecentrality(degree_cent_sf,130);

closs_cent_sf = closeness(sf);
[~,closs_cent_sf1]  = cumulativecentrality(closs_cent_sf,130);

betw_cent_sf  = node_betweenness_faster(sf);
[~,betw_cent_sf1]   = cumulativecentrality(betw_cent_sf,130);

eigen_cent_sf = eigencentrality(full(sf));
[~,eigen_cent_sf1]  = cumulativecentrality(eigen_cent_sf,130);

figure(49);
hold on;
subplot(2,1,1); cdfplot(degree_cent_sf1); title('Deegree centrality SF');
subplot(2,1,2); cdfplot(betw_cent_sf1);   title('Closeness centrality SF');
hold off;
cd StepZ/SF 
print -djpeg CentralitySF.jpeg;
cd ../..

figure(50);
hold on;
subplot(2,1,1); cdfplot(closs_cent_sf1);  title('Betwwness centrality SF');
subplot(2,1,2); cdfplot(eigen_cent_sf1);  title('Eigenvector centality SF');
hold off;
cd StepZ/SF 
print -djpeg CentralitySF1.jpeg;
cd ../..

degree_mean_cent_sf=mean(degree_cent_sf);
closs_avg_cent_sf=mean(closs_cent_sf);
bwtw_avg_cent_sf=mean(betw_cent_sf);
eigen_avg_cent_sf=mean(eigen_cent_sf);


% SW
[degree_cent_sw,~,~] = degrees(full(sw));
[~,degree_cent_sw1]  = cumulativecentrality(degree_cent_sw,130);

closs_cent_sw = closeness(sw);
[~,closs_cent_sw1]  = cumulativecentrality(closs_cent_sw,130);

betw_cent_sw  = node_betweenness_faster(sw);
[~,betw_cent_sw1]   = cumulativecentrality(betw_cent_sw,130);

eigen_cent_sw = eigencentrality(full(sw));
[~,eigen_cent_sw1]  = cumulativecentrality(eigen_cent_sw,130);

figure(51);
hold on;
subplot(2,1,1); cdfplot(degree_cent_sw1); title('Deegree centrality SW');
subplot(2,1,2); cdfplot(betw_cent_sw1);   title('Closeness centrality SW');
hold off;
cd StepZ/SW 
print -djpeg CentralitySW.jpeg;
cd ../..

figure(52);
hold on;
subplot(2,1,1); cdfplot(closs_cent_sw1);  title('Betwwness centrality SW');
subplot(2,1,2); cdfplot(eigen_cent_sw1);  title('Eigenvector centality SW');
hold off;
cd StepZ/SW 
print -djpeg CentralitySW1.jpeg;
cd ../..

degree_mean_cent_sw=mean(degree_cent_sw);
closs_avg_cent_sw=mean(closs_cent_sw);
bwtw_avg_cent_sw=mean(betw_cent_sw);
eigen_avg_cent_sw=mean(eigen_cent_sw);


%% Step H 

mkdir StepH
cd StepH; 
mkdir REG; mkdir ERRG; mkdir RG; mkdir SF; mkdir SW; mkdir RGG;
cd ..

% REG 
ConnectivityPercentage100 = zeros(5,1);
ConnectivityPercentage200 = zeros(5,1);

for i = 1:5
    for j = 1:100
        er = smallw(100,i,0); % to *2 ginetai mesa sthn smallw
        if isconnected(er)
            ConnectivityPercentage100(i) = ConnectivityPercentage100(i)+1;
        end
        er2 = smallw(200,i,0);
        if isconnected(er2)
            ConnectivityPercentage200(i) = ConnectivityPercentage200(i)+1;
        end
    end
    ConnectivityPercentage100(i) = ConnectivityPercentage100(i)/100;
    ConnectivityPercentage200(i) = ConnectivityPercentage200(i)/100;
end

%ConnectivityPercentage100 = ConnectivityPercentage100/100;
%ConnectivityPercentage200 = ConnectivityPercentage200/100;

figure(53);
plot([2:2:10],ConnectivityPercentage100,'blue');
hold on;
plot([2:2:10],ConnectivityPercentage200,'red');
title('Connectivity Percentage of REG to d');
xlabel('Values of d'); ylabel('Connectivity Percentage');
legend('for 100 nodes', 'for 200 nodes','Location', 'SouthEast');
cd StepH/REG 
print -djpeg percentage.jpeg;
cd ../..
hold off;


% ER-RG 
clearvars ConnectivityPercentage100 ConnectivityPercentage200;

ConnectivityPercentage100 = zeros(8,1);
ConnectivityPercentage200 = zeros(8,1);

for i = 1:8
    for j = 1:100
        erd = erdrey(100,i*100);
        if isconnected(erd)
            ConnectivityPercentage100(i) = ConnectivityPercentage100(i)+1;
        end
        erd2 = erdrey(200,i*100);
        if isconnected(erd2)
            ConnectivityPercentage200(i) = ConnectivityPercentage200(i)+1;
        end
    end
end

ConnectivityPercentage100 = ConnectivityPercentage100/100;
ConnectivityPercentage200 = ConnectivityPercentage200/100;

figure(54);
plot([100:100:800],ConnectivityPercentage100,'blue');
hold on;
plot([100:100:800],ConnectivityPercentage200,'red');
title('Connectivity Percentage of ER-RG to M');
xlabel('Values of M'); ylabel('Connectivity Percentage');
legend('for 100 nodes', 'for 200 nodes','Location', 'SouthEast');
cd StepH/ERRG 
print -djpeg percentage.jpeg;
cd ../..
hold off;



%RG
clearvars ConnectivityPercentage100 ConnectivityPercentage200;

ConnectivityPercentage100 = zeros(9,1);
ConnectivityPercentage200 = zeros(9,1);

for i = 1:9
    for j = 1:100
        gil = random_graph(100,0.1*i);
        if isconnected(gil)
            ConnectivityPercentage100(i) = ConnectivityPercentage100(i)+1;
        end
        gil2 = random_graph(200,0.1*i);
        if isconnected(gil2)
            ConnectivityPercentage200(i) = ConnectivityPercentage200(i)+1;
        end
    end
    ConnectivityPercentage100(i) = ConnectivityPercentage100(i)/100;
    ConnectivityPercentage200(i) = ConnectivityPercentage200(i)/100;
end

%ConnectivityPercentage100 = ConnectivityPercentage100/100;
%ConnectivityPercentage200 = ConnectivityPercentage200/100;

figure(55);
plot([0.1:0.1:0.9],ConnectivityPercentage100,'blue');
hold on;
plot([0.1:0.1:0.9],ConnectivityPercentage200,'red');
title('Connectivity Percentage of RG to p');
xlabel('Values of p'); ylabel('Connectivity Percentage');
legend('for 100 nodes', 'for 200 nodes','Location', 'SouthEast');
cd StepH/RG 
print -djpeg percentage.jpeg;
cd ../..
hold off;


%RGG
clearvars ConnectivityPercentage100 ConnectivityPercentage200;

clearvars xrgg yrgg;

for i=1:100
    xrgg(i) = rand()*1000;
    yrgg(i) = rand()*1000;
end
rggnodes100 = [xrgg ; yrgg]';

clearvars xrgg yrgg;

for i=1:200
    xrgg(i) = rand()*1000;
    yrgg(i) = rand()*1000;
end
rggnodes200 = [xrgg ; yrgg]';

ConnectivityPercentage100 = zeros(10,1);
ConnectivityPercentage200 = zeros(10,1);

for i = 1:10
    for j = 1:100
        gil = rgg(rggnodes100,100,i*25);
        if isconnected(gil)
            ConnectivityPercentage100(i) = ConnectivityPercentage100(i)+1;
        end
        gil2 = rgg(rggnodes200,200,i*25);
        if isconnected(gil2)
            ConnectivityPercentage200(i) = ConnectivityPercentage200(i)+1;
        end
    end
end

ConnectivityPercentage100 = ConnectivityPercentage100/100;
ConnectivityPercentage200 = ConnectivityPercentage200/100;

figure(56);
plot([25:25:250],ConnectivityPercentage100,'blue');
hold on;
plot([25:25:250],ConnectivityPercentage200,'red');
title('Connectivity Percentage of RGG to R');
xlabel('Values of R'); ylabel('Connectivity Percentage');
legend('for 100 nodes', 'for 200 nodes','Location', 'SouthEast');
cd StepH/RGG
print -djpeg percentage.jpeg;
cd ../..
hold off;


%BA-SF
clearvars ConnectivityPercentage100 ConnectivityPercentage200;

ConnectivityPercentage100 = zeros(5,1);
ConnectivityPercentage200 = zeros(5,1);

for i = 1:5
    for j = 1:100
        tmp = pref(100,i*2);;
        if isconnected(tmp)
            ConnectivityPercentage100(i) = ConnectivityPercentage100(i)+1;
        end
        tmp2 = pref(200,i*2);;
        if isconnected(tmp2)
            ConnectivityPercentage200(i) = ConnectivityPercentage200(i)+1;
        end
    end
    ConnectivityPercentage100(i) = ConnectivityPercentage100(i)/100;
    ConnectivityPercentage200(i) = ConnectivityPercentage200(i)/100;
end

figure(57);
plot([2:2:10],ConnectivityPercentage100,'blue');
hold on;
plot([2:2:10],ConnectivityPercentage200,'red');
title('Connectivity Percentage of BA-SF to d');
xlabel('Values of d'); ylabel('Connectivity Percentage');
legend('for 100 nodes', 'for 200 nodes','Location', 'SouthEast');
cd StepH/SF 
print -djpeg percentage.jpeg;
cd ../..
hold off;

%SW 
clearvars ConnectivityPercentage100 ConnectivityPercentage200;

ConnectivityPercentage100 = zeros(5,1);
ConnectivityPercentage200 = zeros(5,1);

for i = 1:5
    for j = 1:100
        er = smallw(100,i,0.1); % to *2 ginetai mesa sthn smallw
        if isconnected(er)
            ConnectivityPercentage100(i) = ConnectivityPercentage100(i)+1;
        end
        er2 = smallw(200,i,0.1);
        if isconnected(er2)
            ConnectivityPercentage200(i) = ConnectivityPercentage200(i)+1;
        end
    end
    ConnectivityPercentage100(i) = ConnectivityPercentage100(i)/100;
    ConnectivityPercentage200(i) = ConnectivityPercentage200(i)/100;
end

%ConnectivityPercentage100 = ConnectivityPercentage100/100;
%ConnectivityPercentage200 = ConnectivityPercentage200/100;

figure(58);
plot([2:2:10],ConnectivityPercentage100,'blue');
hold on;
plot([2:2:10],ConnectivityPercentage200,'red');
title('Connectivity Percentage of SW to d with g=0.1');
xlabel('Values of d'); ylabel('Connectivity Percentage');
legend('for 100 nodes', 'for 200 nodes','Location', 'SouthEast');
cd StepH/SW 
print -djpeg percentage.jpeg;
cd ../..
hold off;

clearvars ConnectivityPercentage100 ConnectivityPercentage200;

ConnectivityPercentage100 = zeros(7,1);
ConnectivityPercentage200 = zeros(7,1);

for i = 1:7
    for j = 1:100
        er = smallw(100,2,0.1*i); % to *2 ginetai mesa sthn smallw
        if isconnected(er)
            ConnectivityPercentage100(i) = ConnectivityPercentage100(i)+1;
        end
        er2 = smallw(200,2,0.1*i);
        if isconnected(er2)
            ConnectivityPercentage200(i) = ConnectivityPercentage200(i)+1;
        end
    end
    ConnectivityPercentage100(i) = ConnectivityPercentage100(i)/100;
    ConnectivityPercentage200(i) = ConnectivityPercentage200(i)/100;
end

ConnectivityPercentage100 = ConnectivityPercentage100/100;
ConnectivityPercentage200 = ConnectivityPercentage200/100;

figure(59);
plot([0.1:0.1:0.7],ConnectivityPercentage100,'blue');
hold on;
plot([0.1:0.1:0.7],ConnectivityPercentage200,'red');
title('Connectivity Percentage of SW to g with d=2');
xlabel('Values of g'); ylabel('Connectivity Percentage');
legend('for 100 nodes', 'for 200 nodes','Location', 'SouthEast');
cd StepH/SW 
print -djpeg percentage.jpeg;
cd ../..
hold off;

%% Step I 

j=1;
for i=0:0.1:1
    sw_I=smallw(130,2,i);
    ave_path_sw_I(j) = ave_path_length(sw_I);
    [~,average_clust_sw_I(j),~] = clust_coeff(sw_I);
    j=j+1;
end

%% Step K - Ego Centralities
mkdir StepK
cd StepK; 
mkdir REG; mkdir RG; mkdir SF; mkdir SW; mkdir RGG;
cd ..

% REG
ego_centr = zeros(130,1);
for i=1:130
    T = [i,kneighbors(regA,i,1)]; 
    A = subgraph(regA,T);
    ego_A = A*A*(ones(size(A,1),size(A,2))-A);
    ego_A = triu(ego_A,1);
    ego_A = ego_A .* (1-A);
    ego_A = ego_A(ego_A~=0);
    ego_A = 1./ego_A;
    ego_centr(i) = sum(sum(ego_A));
end 

[ego_sort,I_e] = sort(ego_centr);
[bet_sort,I_b] = sort(betw_cent_reg);

result = zeros(130,1);
for i = 1 : 130
    ego_point = find(I_e == i);
    bet_point = find(I_b == i);
    result(i) = abs(ego_point - bet_point);
end

figure(60);
stem(ego_sort,'filled','b');  
set(gca, 'XTick', 1:length(I_e)); % Change x-axis ticks
set(gca, 'XTickLabel', I_e); % Change x-axis ticks labels.
ax = gca; ax.XTickLabelRotation=90;
cd StepK/REG 
title('Sorted Nodes of REG Ego Betweenness');
print -djpeg egonodes.jpeg;
cd ../..

figure(61);
stem(bet_sort,'filled','b');  
set(gca, 'XTick', 1:length(I_b)); % Change x-axis ticks
set(gca, 'XTickLabel', I_b); % Change x-axis ticks labels.
ax = gca; ax.XTickLabelRotation=90;
cd StepK/REG 
title('Sorted Nodes of REG Betweenness');
print -djpeg betnodes.jpeg;
cd ../..

figure(62);
stem(result,'filled','b');
ylabel('Abs Distace');
xlabel('Nodes from 1 to 130')
cd StepK/REG 
title('Abs of nodes positon distances');
print -djpeg nodesdistance.jpeg;
cd ../..


% RG
clearvars ego_centr T A ego_A ego_sort I_e I_b bet_sort

ego_centr = zeros(130,1);
for i=1:130
    T = [i,kneighbors(rgA,i,1)]; 
    A = subgraph(rgA,T);
    ego_A = A*A*(ones(size(A,1),size(A,2))-A);
    ego_A = triu(ego_A,1);
    ego_A = ego_A .* (1-A);
    ego_A = ego_A(ego_A~=0);
    ego_A = 1./ego_A;
    ego_centr(i) = sum(sum(ego_A));
end 

[ego_sort,I_e] = sort(ego_centr);
[bet_sort,I_b] = sort(betw_cent_rg);

result = zeros(130,1);
for i = 1 : 130
    ego_point = find(I_e == i);
    bet_point = find(I_b == i);
    result(i) = abs(ego_point - bet_point);
end

figure(63);
stem(ego_sort,'filled','b');  
set(gca, 'XTick', 1:length(I_e)); % Change x-axis ticks
set(gca, 'XTickLabel', I_e); % Change x-axis ticks labels.
ax = gca; ax.XTickLabelRotation=90;
cd StepK/RG 
title('Sorted Nodes of RG Ego Betweenness');
print -djpeg egonodes.jpeg;
cd ../..

figure(64);
stem(bet_sort,'filled','b');  
set(gca, 'XTick', 1:length(I_b)); % Change x-axis ticks
set(gca, 'XTickLabel', I_b); % Change x-axis ticks labels.
ax = gca; ax.XTickLabelRotation=90;
cd StepK/RG 
title('Sorted Nodes of RG Betweenness');
print -djpeg betnodes.jpeg;
cd ../..

figure(65);
stem(result,'filled','b');
ylabel('Abs Distace');
xlabel('Nodes from 1 to 130')
cd StepK/RG 
title('Abs of nodes positon distances');
print -djpeg nodesdistance.jpeg;
cd ../..


% RGG
clearvars ego_centr T A ego_A ego_sort I_e I_b bet_sort

ego_centr = zeros(130,1);
for i=1:130
    T = [i,kneighbors(rggA,i,1)]; 
    A = subgraph(rggA,T);
    ego_A = A*A*(ones(size(A,1),size(A,2))-A);
    ego_A = triu(ego_A,1);
    ego_A = ego_A .* (1-A);
    ego_A = ego_A(ego_A~=0);
    ego_A = 1./ego_A;
    ego_centr(i) = sum(sum(ego_A));
end 

[ego_sort,I_e] = sort(ego_centr);
[bet_sort,I_b] = sort(betw_cent_rgg);

result = zeros(130,1);
for i = 1 : 130
    ego_point = find(I_e == i);
    bet_point = find(I_b == i);
    result(i) = abs(ego_point - bet_point);
end

figure(66);
stem(ego_sort,'filled','b');  
set(gca, 'XTick', 1:length(I_e)); % Change x-axis ticks
set(gca, 'XTickLabel', I_e); % Change x-axis ticks labels.
ax = gca; ax.XTickLabelRotation=90;
cd StepK/RGG
title('Sorted Nodes of RG Ego Betweenness');
print -djpeg egonodes.jpeg;
cd ../..

figure(67);
stem(bet_sort,'filled','b');  
set(gca, 'XTick', 1:length(I_b)); % Change x-axis ticks
set(gca, 'XTickLabel', I_b); % Change x-axis ticks labels.
ax = gca; ax.XTickLabelRotation=90;
cd StepK/RGG
title('Sorted Nodes of RGG Betweenness');
print -djpeg betnodes.jpeg;
cd ../..

figure(68);
stem(result,'filled','b');
ylabel('Abs Distace');
xlabel('Nodes from 1 to 130')
cd StepK/RGG
title('Abs of nodes positon distances');
print -djpeg nodesdistance.jpeg;
cd ../..


% SF
clearvars ego_centr T A ego_A ego_sort I_e I_b bet_sort

ego_centr = zeros(130,1);
for i=1:130
    T = [i,kneighbors(sf,i,1)]; 
    A = subgraph(sf,T);
    ego_A = A*A*(ones(size(A,1),size(A,2))-A);
    ego_A = triu(ego_A,1);
    ego_A = ego_A .* (1-A);
    ego_A = ego_A(ego_A~=0);
    ego_A = 1./ego_A;
    ego_centr(i) = sum(sum(ego_A));
end 

[ego_sort,I_e] = sort(ego_centr);
[bet_sort,I_b] = sort(betw_cent_sf);

result = zeros(130,1);
for i = 1 : 130
    ego_point = find(I_e == i);
    bet_point = find(I_b == i);
    result(i) = abs(ego_point - bet_point);
end

figure(69);
stem(ego_sort,'filled','b');  
set(gca, 'XTick', 1:length(I_e)); % Change x-axis ticks
set(gca, 'XTickLabel', I_e); % Change x-axis ticks labels.
ax = gca; ax.XTickLabelRotation=90;
cd StepK/SF
title('Sorted Nodes of SF Ego Betweenness');
print -djpeg egonodes.jpeg;
cd ../..

figure(70);
stem(bet_sort,'filled','b');  
set(gca, 'XTick', 1:length(I_b)); % Change x-axis ticks
set(gca, 'XTickLabel', I_b); % Change x-axis ticks labels.
ax = gca; ax.XTickLabelRotation=90;
cd StepK/SF
title('Sorted Nodes of SF Betweenness');
print -djpeg betnodes.jpeg;
cd ../..

figure(71);
stem(result,'filled','b');
ylabel('Abs Distace');
xlabel('Nodes from 1 to 130')
cd StepK/SF
title('Abs of nodes positon distances');
print -djpeg nodesdistance.jpeg;
cd ../..


% SW
clearvars ego_centr T A ego_A ego_sort I_e I_b bet_sort

ego_centr = zeros(130,1);
for i=1:130
    T = [i,kneighbors(sw,i,1)]; 
    A = subgraph(sw,T);
    ego_A = A*A*(ones(size(A,1),size(A,2))-A);
    ego_A = triu(ego_A,1);
    ego_A = ego_A .* (1-A);
    ego_A = ego_A(ego_A~=0);
    ego_A = 1./ego_A;
    ego_centr(i) = sum(sum(ego_A));
end 

[ego_sort,I_e] = sort(ego_centr);
[bet_sort,I_b] = sort(betw_cent_sw);

result = zeros(130,1);
for i = 1 : 130
    ego_point = find(I_e == i);
    bet_point = find(I_b == i);
    result(i) = abs(ego_point - bet_point);
end

figure(72);
stem(ego_sort,'filled','b');  
set(gca, 'XTick', 1:length(I_e)); % Change x-axis ticks
set(gca, 'XTickLabel', I_e); % Change x-axis ticks labels.
ax = gca; ax.XTickLabelRotation=90;
cd StepK/SW
title('Sorted Nodes of SW Ego Betweenness');
print -djpeg egonodes.jpeg;
cd ../..

figure(73);
stem(bet_sort,'filled','b');  
set(gca, 'XTick', 1:length(I_b)); % Change x-axis ticks
set(gca, 'XTickLabel', I_b); % Change x-axis ticks labels.
ax = gca; ax.XTickLabelRotation=90;
cd StepK/SW
title('Sorted Nodes of SW Betweenness');
print -djpeg betnodes.jpeg;
cd ../..

figure(74);
stem(result,'filled','b');
ylabel('Abs Distace');
xlabel('Nodes from 1 to 130')
cd StepK/SW
title('Abs of nodes positon distances');
print -djpeg nodesdistance.jpeg;
cd ../..
