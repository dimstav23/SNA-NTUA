%% SNA Lab2 
clear all; close all; clc;
PathAdd
tic;
%% Step A 
mkdir StepA;
cd StepA; mkdir football; mkdir dolphins; mkdir lesmis;
cd ..
disp('************ <Step A> ************');

lesmis   = importgml('lesmis.gml');
dolphins = importgml('dolphins.gml');
football = importgml('football.gml');

% check if graphs are directed or not and if it is make it undirected
if isdirected(football)
    temp=undir(football,size(football,1));
    football=temp;
    clearvars temp;
end
if isdirected(lesmis)
    temp=undir(football,size(lesmis,1));
    lesmis=temp;
    clearvars temp;
end
if isdirected(dolphins)
    temp=undir(football,size(dolphins,1));
    dolphins=temp;
    clearvars temp;
end

% Football Results
disp('Football Results');
foot_deg = degrees(football); % vathmos ka8e komvou
foot_mean_deg = mean(foot_deg);
disp(['football mean degree:', num2str(foot_mean_deg)]);

ego_foot = MyEgoCent(football,115); % EGO
mean_ego_foot = mean(ego_foot);
ego_foot_2 = ego_cent(football,115);
mean_ego_foot_2 = mean(ego_foot_2);
ego_foot_3 = new_ego(football,115);
mean_ego_foot_3 = mean(ego_foot_3);

disp(['football mean ego:', num2str(mean_ego_foot)]);

[~,foot_c_mean,foot_clust] = clust_coeff(full(football)); % Clust Coeff 
disp(['football mean cc:',num2str(foot_c_mean)]);

figure(1);
set(gcf, 'Position', get(0,'Screensize'));
a = subplot(3,1,3);
stem(1:115,ego_foot,'MarkerFaceColor','r','Color','c');
title(a,'Ego Centrality');
b = subplot(3,1,1);
stem(1:115,foot_deg','MarkerFaceColor','r','Color','g');
title(b,'Degree of each node');
c = subplot(3,1,2);
stem(1:115,foot_clust,'MarkerFaceColor','r','Color','y');
title(c,'Clustering Coefficient');
suptitle('Football');
cd StepA/football
print -djpeg Football.jpeg;
cd ../..

% Lesmis Results
disp('Lesmis Results');
les_deg = degrees(lesmis); % vathmos ka8e komvou
les_mean_deg = mean(les_deg);
disp(['lesmis   mean degree:', num2str(les_mean_deg)]);

ego_les = MyEgoCent(lesmis,77); % EGO
mean_ego_les = mean(ego_les);
ego_les_2 = ego_cent(lesmis,77);
mean_ego_les_2 = mean(ego_les_2);
disp(['lesmis   mean ego:', num2str(mean_ego_les)]);

[~,les_c_mean,les_clust] = clust_coeff(full(lesmis)); % Clust Coeff 
disp(['lesmis   mean cc:',num2str(les_c_mean)]);

figure(2);
set(gcf, 'Position', get(0,'Screensize'));
a = subplot(3,1,3);
stem(1:77,ego_les,'MarkerFaceColor','r','Color','c');
title(a,'Ego Centrality');
b = subplot(3,1,1);
stem(1:77,les_deg','MarkerFaceColor','r','Color','g');
title(b,'Degree of each node');
c = subplot(3,1,2);
stem(1:77,les_clust,'MarkerFaceColor','r','Color','y');
title(c,'Clustering Coefficient');
suptitle('Les Miserables');
cd StepA/lesmis
print -djpeg Lesmis.jpeg;
cd ../..

% Dolphins Results
disp('Dolphins Results');
dol_deg = degrees(dolphins); % vathmos ka8e komvou
dol_mean_deg = mean(dol_deg);
disp(['dolphins mean degree:', num2str(dol_mean_deg)]);

ego_dol = MyEgoCent(dolphins,62); % EGO
mean_ego_dol = mean(ego_dol);
ego_dol_2 = ego_cent(dolphins,62);
mean_ego_dol_2 = mean(ego_dol_2);
disp(['dolphins mean ego:', num2str(mean_ego_dol)]);

[~,dol_c_mean,dol_clust] = clust_coeff(full(dolphins)); % Clust Coeff 
disp(['dolphins mean cc:',num2str(dol_c_mean)]);

figure(3);
set(gcf, 'Position', get(0,'Screensize'));
a = subplot(3,1,3);
stem(1:62,ego_dol,'MarkerFaceColor','r','Color','c');
title(a,'Ego Centrality');
b = subplot(3,1,1);
stem(1:62,dol_deg','MarkerFaceColor','r','Color','g');
title(b,'Degree of each node');
c = subplot(3,1,2);
stem(1:62,dol_clust,'MarkerFaceColor','r','Color','y');
title(c,'Clustering Coefficient');
suptitle('Dolphins');
cd StepA/dolphins
print -djpeg Dolphins.jpeg;
cd ../..

disp('************ <End of Step A> ************');
close all;
%% Step B  
mkdir StepB;
cd StepB; mkdir football; mkdir dolphins; mkdir lesmis; mkdir REG; mkdir RG; mkdir RGG; mkdir SF; mkdir SW;
cd ..
disp('************ <Step B> ************');

% REG
regA=smallw(170,2,0);
[newman_reg, newman_modul_reg]   = newman_comm_fast(regA);
[newman_max_reg,newman_max_reg1] = max(newman_modul_reg);

figure(4);
PlotGraph(full(regA),clustercell(newman_reg,170,newman_max_reg1));
title('REG:Newman-Girvan method');
cd StepB/REG
print -djpeg newman.jpeg;
cd ../..
reg_spectral = GCSpectralClust2(regA,170,4);
spectralQ = zeros(170,1);
for i=1:170
    spectralQ(i)=QFModul(reg_spectral(:,i),regA);
end
[spectral_max_reg,spectral_max_reg1]=max(spectralQ);
figure(5);
PlotGraph(full(regA),reg_spectral(:,spectral_max_reg1));
title('REG:Spectral Clustering method');
cd StepB/REG
print -djpeg spectral.jpeg;
cd ../..
Modulemax_reg = GCModulMax1(regA);
modulemax_val_reg = QFModul(Modulemax_reg,regA);
figure(6);
PlotGraph(full(regA),Modulemax_reg);
title('REG:Modularity Maximization method');
cd StepB/REG
print -djpeg maxmodule.jpeg;
cd ../..

% RG
clearvars SpeactralQ;
rgA=erdrey(170,750);
[newman_rg, newman_modul_rg]= newman_comm_fast(rgA);
[newman_max_rg,newman_max_rg1]=max(newman_modul_rg);
figure(6);
PlotGraph(full(rgA),clustercell(newman_rg,170,newman_max_rg1));
title('RG:Newman-Girvan method');
cd StepB/RG
saveas(gcf,'Newman_rg.png');
cd ../..
rg_spectral=GCSpectralClust2(rgA,170,4);
spectralQ=zeros(170,1);
for i=1:170
    spectralQ(i)=QFModul(rg_spectral(:,i),rgA);
end
[spectral_max_rg,spectral_max_rg1]=max(spectralQ);
figure(7);
PlotGraph(full(rgA),rg_spectral(:,spectral_max_rg1));
title('RG:Spectral Clustering method');
cd StepB/RG
saveas(gcf,'Spectral_RG.png');
cd ../..
Modulemax_rg=GCModulMax1(rgA);
modulemax_val_rg= QFModul(Modulemax_rg,rgA);
figure(8);
PlotGraph(full(rgA),Modulemax_rg);
title('RG:Modularity Maximization method');
cd StepB/RG
saveas(gcf,'Module_rg.png');
cd ../..

%% RGG
clearvars x y;
%generates random points in LxL
for i=1:180
    x(i)=rand()*1000; % a random real number between 0-1 * norma . like a pointer
    y(i)=rand()*1000;
end
rgNodes = transpose([x ; y]);
[rggA,rggNodeDegree]=rgg(rgNodes,170,250);
clearvars SpeactralQ;
[newman_rgg, newman_modul_rgg]= newman_comm_fast(rggA);
[newman_max_rgg,newman_max_rgg1]=max(newman_modul_rgg);
figure(9);
PlotGraph(full(rggA),clustercell(newman_rgg,170,newman_max_rgg1));
title('RGG:Newman-Girvan method');
cd StepB/RGG
saveas(gcf,'Newman_rgg.png');
cd ../..
rgg_spectral=GCSpectralClust2(rggA,170,8);
spectralQ=zeros(170,1);
for i=1:170
    spectralQ(i)=QFModul(rgg_spectral(:,i),rggA);
end
[spectral_max_rgg,spectral_max_rgg1]=max(spectralQ);
figure(10);
PlotGraph(full(rggA),rgg_spectral(:,spectral_max_rgg1));
title('RGG:Spectral Clustering method');
cd StepB/RGG
saveas(gcf,'Spectral_RGG.png');
cd ../..
Modulemax_rgg=GCModulMax1(rggA);
modulemax_val_rgg= QFModul(Modulemax_rgg,rggA);
figure(11);
PlotGraph(full(rggA),Modulemax_rgg);
title('RGG:Modularity Maximization method');
cd StepB/RGG
saveas(gcf,'Module_rgg.png');
cd ../..

%% SF
sf=pref(170,4);
clearvars SpeactralQ;
[newman_sf, newman_modul_sf]= newman_comm_fast(sf);
[newman_max_sf,newman_max_sf1]=max(newman_modul_sf);
figure(12);
PlotGraph(full(sf),clustercell(newman_sf,170,newman_max_sf1));
title('SF:Newman-Girvan method');
cd StepB/SF
saveas(gcf,'Newman_sf.png');
cd ../..
sf_spectral=GCSpectralClust2(sf,170,4);
spectralQ=zeros(170,1);
for i=1:170
    spectralQ(i)=QFModul(sf_spectral(:,i),sf);
end
[spectral_max_sf,spectral_max_sf1]=max(spectralQ);
figure(13);
PlotGraph(full(sf),sf_spectral(:,spectral_max_sf1));
title('SF:Spectral Clustering method');
cd StepB/SF
saveas(gcf,'Spectral_SF.png');
cd ../..
Modulemax_sf=GCModulMax1(sf);
modulemax_val_sf= QFModul(Modulemax_sf,sf);
figure(14);
PlotGraph(full(sf),Modulemax_sf);
title('SF:Modularity Maximization method');
cd StepB/SF
saveas(gcf,'Module_sf.png');
cd ../..

%% SW
sw=smallw(170,2,0.3);
%foot_deg = degrees(sw); % vathmos ka8e komvou
%stem(1:170,foot_deg','MarkerFaceColor','r','Color','y');
clearvars SpeactralQ;
[newman_sw, newman_modul_sw]= newman_comm_fast(sw);
[newman_max_sw,newman_max_sw1]=max(newman_modul_sw);
figure(15);
PlotGraph(full(sw),clustercell(newman_sw,180,newman_max_sw1));
title('SW:Newman-Girvan method');
cd StepB/SW
saveas(gcf,'Newman_sw.png');
cd ../..
sw_spectral=GCSpectralClust2(sw,170,4);
spectralQ=zeros(170,1);
for i=1:170
    spectralQ(i)=QFModul(sw_spectral(:,i),sw);
end
[spectral_max_sw,spectral_max_sw1]=max(spectralQ);
figure(16);
PlotGraph(full(sw),sw_spectral(:,spectral_max_sw1));
title('SW:Spectral Clustering method');
cd StepB/SW
saveas(gcf,'Spectral_SW.png');
cd ../..
Modulemax_sw=GCModulMax1(sw);
modulemax_val_sw= QFModul(Modulemax_sw,sw);
figure(17);
PlotGraph(full(sw),Modulemax_sw);
title('SW:Modularity Maximization method');
cd StepB/SW
saveas(gcf,'Module_sw.png');
cd ../..

% Dolphins

clearvars SpeactralQ;
[newman_dolphins, newman_modul_dolphins]= newman_comm_fast(dolphins);
[newman_max_dolphins,newman_max_dolphins1]=max(newman_modul_dolphins);
figure(18);
PlotGraph(full(dolphins),clustercell(newman_dolphins,62,newman_max_dolphins1));
title('Dolphins Graph:Newman-Girvan method');
cd StepB/dolphins
saveas(gcf,'Newman_dolphins.png');
cd ../..
dolphins_spectral=GCSpectralClust2(dolphins,62,4);
spectralQ=zeros(62,1);
for i=1:62
    spectralQ(i)=QFModul(dolphins_spectral(:,i),dolphins);
end
[spectral_max_dolphins,spectral_max_dolphins1]=max(spectralQ);
figure(19);
PlotGraph(full(dolphins),dolphins_spectral(:,spectral_max_dolphins1));
title('Dolphins Graph:Spectral Clustering method');
cd StepB/dolphins
saveas(gcf,'Spectral_dolphins.png');
cd ../..
Modulemax_dolphins=GCModulMax1(dolphins);
modulemax_val_dolphins= QFModul(Modulemax_dolphins,dolphins);
figure(20);
PlotGraph(full(dolphins),Modulemax_dolphins);
title('Dolphins Graph:Modularity Maximization method');
cd StepB/dolphins
saveas(gcf,'Module_dolphins.png');
cd ../..

% lesmis
clearvars SpeactralQ;
[newman_lesmis, newman_modul_lesmis]= newman_comm_fast(lesmis);
[newman_max_lesmis,newman_max_lesmis1]=max(newman_modul_lesmis);
figure(21);
PlotGraph(full(lesmis),clustercell(newman_lesmis,77,newman_max_lesmis1));
title('Les Miserables Graph:Newman-Girvan method');
cd StepB/lesmis
saveas(gcf,'Newman_lesmis.png');
cd ../..
lesmis_spectral=GCSpectralClust2(lesmis,77,4);
spectralQ=zeros(77,1);
for i=1:77
    spectralQ(i)=QFModul(lesmis_spectral(:,i),lesmis);
end
[spectral_max_lesmis,spectral_max_lesmis1]=max(spectralQ);
figure(22);
PlotGraph(full(lesmis),lesmis_spectral(:,spectral_max_lesmis1));
title('Les Miserables Graph:Spectral Clustering method');
cd StepB/lesmis
saveas(gcf,'Spectral_lesmis.png');
cd ../..
Modulemax_lesmis=GCModulMax1(lesmis);
modulemax_val_lesmis= QFModul(Modulemax_lesmis,lesmis);
figure(22);
PlotGraph(full(lesmis),Modulemax_lesmis);
title('Les Miserables Graph:Modularity Maximization method');
cd StepB/lesmis
saveas(gcf,'Module_lesmis.png');
cd ../..


% Football
clearvars SpeactralQ;
[newman_football, newman_modul_football]= newman_comm_fast(football);
[newman_max_football,newman_max_football1]=max(newman_modul_football);
figure(23);
PlotGraph(full(football),clustercell(newman_football,62,newman_max_football1));
title('Football Graph:Newman-Girvan method');
cd StepB/football
saveas(gcf,'Newman_football.png');
cd ../..
football_spectral=GCSpectralClust2(football,62,4);
spectralQ=zeros(62,1);
for i=1:62
    spectralQ(i)=QFModul(football_spectral(:,i),football);
end
[spectral_max_football,spectral_max_football1]=max(spectralQ);
figure(24);
PlotGraph(full(football),football_spectral(:,spectral_max_football1));
title('Football Graph:Spectral Clustering method');
cd StepB/football
saveas(gcf,'Spectral_football.png');
cd ../..
Modulemax_football=GCModulMax1(football);
modulemax_val_football= QFModul(Modulemax_football,football);
figure(25);
PlotGraph(full(football),Modulemax_football);
title('Football Graph:Modularity Maximization method');
cd StepB/football
saveas(gcf,'Module_football.png');
cd ../..

toc;

