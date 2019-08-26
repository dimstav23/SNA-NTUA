regA=smallw(180,2,0);
[x,y]= getNodeCoordinates(180);
regNodes=[x y];
PathAdd


rgA=erdrey(180,750);

clearvars x y;
%generates random points in LxL
for i=1:180
    x(i)=rand()*1000; % a random real number between 0-1 * norma . like a pointer
    y(i)=rand()*1000;
end

rgNodes = transpose([x ; y]);

[rggA,rggNodeDegree]=rgg(rgNodes,180,250);



sf=pref(180,4);


sw=smallw(180,2,0.3); % the difference with the first is that now we have set the propability.
%% START OF EXERISE 2
%% EGO  REG
for i=1:180
    T=[i,kneighbors(regA,i,1)];
    temp= subgraph(regA,T);
    len=size(full(T));
    ego_reg(i)= ego_net(full(temp),len(2));
    clearvars len temp T;
end
figure ;
s = subplot(1,1,1);
stem(s,[1:180],ego_reg,'filled','r');
title('EGO REG');
bet_reg=brandesBetwCentr(regA)';
reg_cmp_r = [ego_reg,1:180]';
reg_cmp_r =sortrows(reg_cmp_r,1);
reg_cmp_d = [bet_reg,1:180]';
reg_cmp_d =sortrows(reg_cmp_d,1);
reg_cmp_t= [reg_cmp_r,reg_cmp_d];

%%
%% EGO  RG


for i=1:180
    T=[i,kneighbors(rgA,i,1)];
    temp= subgraph(rgA,T);
    len=size(full(T));
    ego_rg(i)= ego_net(full(temp),len(2));
    clearvars len temp T;
end
figure ;
s = subplot(1,1,1);
stem(s,[1:180],ego_rg,'filled','r');
title('EGO RG');
saveas(gcf,'ego_rg.png');
bet_rg=brandesBetwCentr(rgA)';
rg_cmp_r = [ego_rg;1:180]';
rg_cmp_r =sortrows(rg_cmp_r,1);
rg_cmp_d = [bet_rg;1:180]';
rg_cmp_d =sortrows(rg_cmp_d,1);
rg_cmp_t = [rg_cmp_r , rg_cmp_d];
%% EGO rgg
for i=1:180
    T=[i,kneighbors(rggA,i,1)];
    temp= subgraph(rggA,T);
    len=size(full(T));
    ego_rgg(i)= ego_net(full(temp),len(2));
    clearvars len temp T;
end
figure ;
s = subplot(1,1,1);
stem(s,[1:180],ego_rgg,'filled','r');
title('EGO RGG');
saveas(gcf,'ego_rgg.png');
bet_rgg=brandesBetwCentr(rggA)';
rgg_cmp_r = [ego_rgg;1:180]';
rgg_cmp_r =sortrows(rgg_cmp_r,1);
rgg_cmp_d = [bet_rgg;1:180]';
rgg_cmp_d =sortrows(rgg_cmp_d,1);
rgg_cmp_t= [rgg_cmp_r,rgg_cmp_d];
%% EGO SF
for i=1:180
    T=[i,kneighbors(full(sf),i,1)];
    temp= subgraph(full(sf),T);
    len=size(full(T));
    ego_sf(i)= ego_net(full(temp),len(2));
    clearvars len temp T;
end
figure ;
s = subplot(1,1,1);
stem(s,[1:180],ego_sf,'filled','r');
bet_sf=brandesBetwCentr(full(sf))';
title('EGO SF');
saveas(gcf,'ego_sf.png');
sf_cmp_r = [ego_sf;1:180]';
sf_cmp_r =sortrows(sf_cmp_r,1);
sf_cmp_d = [bet_sf;1:180]';
sf_cmp_d =sortrows(sf_cmp_d,1);
sf_cmp_t= [sf_cmp_r,sf_cmp_d];
%% EGO SW
for i=1:180
    T=[i,kneighbors(sw,i,1)];
    temp= subgraph(sw,T);
    len=size(full(T));
    ego_sw(i)= ego_net(full(temp),len(2));
    clearvars len temp T;
end
figure ;
s = subplot(1,1,1);
stem(s,[1:180],ego_sw,'filled','r');
title('EGO SW');
saveas(gcf,'ego_sw.png');
bet_sw=brandesBetwCentr(sw)';
sw_cmp_r = [ego_sw;1:180]';
sw_cmp_r =sortrows(sw_cmp_r,1);
sw_cmp_d = [bet_sw;1:180]';
sw_cmp_d =sortrows(sw_cmp_d,1);
sw_cmp_t= [sw_cmp_r,sw_cmp_d];
%% END EX1
%% EX2

clearvars x y;
%generates random points in LxL
for i=1:180
    x(i)=rand()*1000; % a random real number between 0-1 * norma . like a pointer
    y(i)=rand()*1000;
end

rgNodes = transpose([x ; y]);

[football,rggNodeDegree]=rgg(rgNodes,115,250);

[lesmis,rggNodeDegree]=rgg(rgNodes,77,250);
[dolphins,rggNodeDegree]=rgg(rgNodes,62,250);
%lesmis =smallw(77,2,0.5); %importgml('lesmis.gml');
%dolphins = smallw(62,2,0.5); %importgml('dolphins.gml');
%Lets make the graphs matrix
if isdirected(football)
    B=undir(football,size(football,1));
    football=B;
    clearvars B;
end
if isdirected(lesmis)
    B=undir(football,size(lesmis,1));
    lesmis=B;
    clearvars B;
end
if isdirected(dolphins)
    B=undir(football,size(dolphins,1));
    dolphins=B;
    clearvars B;
end
ego_foot=ego_cent_all(football,115);
mean_ego_foot= mean(ego_foot);
foot_deg=degrees(football);
foot_mean_deg=mean(foot_deg);
[~,foot_c_mean,foot_clust]= clust_coeff(full(football));
figure;
subplot(3,1,1);
stem(1:115,ego_foot,'filled','g');
subplot(3,1,2);
stem(1:115,foot_deg','filled');
subplot(3,1,3);
stem(1:115,foot_clust,'filled','r');
title('Football: Ego(G),deg(B),clust(R)');
%% dolphins
ego_dolphins=ego_cent_all(dolphins,62);
mean_ego_doplhins= mean(ego_dolphins);
dolphins_deg=degrees(dolphins);
dol_mean_deg=mean(dolphins_deg);
[~,dol_c_mean,dolphins_clust]= clust_coeff(full(dolphins));
figure;
subplot(3,1,1);
stem(1:62,ego_dolphins,'filled','g');
subplot(3,1,2);
stem(1:62,dolphins_deg','filled');
subplot(3,1,3);
stem(1:62,dolphins_clust,'filled','r');
title('Dolphins: Ego(G),deg(B),clust(R)');
%%
% lesmis
ego_lesmis=ego_cent_all(lesmis,77);
mean_ego_lesmis= mean(ego_lesmis);
lesmis_deg=degrees(lesmis);
lesmis_mean_deg=mean(lesmis_deg);
[~,lesmis_c_mean,lesmis_clust]= clust_coeff(full(lesmis));
figure;
subplot(3,1,1);
stem(1:77,ego_lesmis,'filled','g');
subplot(3,1,2);
stem(1:77,lesmis_deg','filled');
subplot(3,1,3);
stem(1:77,lesmis_clust,'filled','r');
title('Lesmis: Ego(G),deg(B),clust(R)');
%% EX3
%% REG
[newman_reg, newman_modul_reg]= newman_comm_fast(regA);
[newman_max_reg,newman_max_reg1]=max(newman_modul_reg);
%figure;
%plot(newman_modul_reg);
%title('Newman Module');
figure;
PlotGraph(full(regA),clustercell(newman_reg,180,newman_max_reg1));
title('Newman method communities REG');
reg_spectral=GCSpectralClust2(regA,180,4);
spectralQ=zeros(180,1);
for i=1:180
    spectralQ(i)=QFModul(reg_spectral(:,i),regA);
end
[spectral_max_reg,spectral_max_reg1]=max(spectralQ);
figure;
PlotGraph(full(regA),reg_spectral(:,spectral_max_reg1));
title('Spectral method communities REG');

Modulemax_reg=GCModulMax1(regA);
modulemax_val_reg= QFModul(Modulemax_reg,regA);
figure;
PlotGraph(full(regA),Modulemax_reg);
title('Max Module method communities REG');

%% RG
clearvars SpeactralQ;
[newman_rg, newman_modul_rg]= newman_comm_fast(rgA);
[newman_max_rg,newman_max_rg1]=max(newman_modul_rg);
%figure;
%plot(newman_modul_reg);
%title('Newman Module');
figure;
PlotGraph(full(rgA),clustercell(newman_rg,180,newman_max_rg1));
title('Newman method communities RG');
saveas(gcf,'Newman_rg.png');
rg_spectral=GCSpectralClust2(rgA,180,4);
spectralQ=zeros(180,1);
for i=1:180
    spectralQ(i)=QFModul(rg_spectral(:,i),rgA);
end
[spectral_max_rg,spectral_max_rg1]=max(spectralQ);
figure;
PlotGraph(full(rgA),rg_spectral(:,spectral_max_rg1));
title('Spectral method communities RG');
saveas(gcf,'Spectral_RG.png');
Modulemax_rg=GCModulMax1(rgA);
modulemax_val_rg= QFModul(Modulemax_rg,rgA);
figure;
PlotGraph(full(rgA),Modulemax_rg);
title('Max Module method communities RG');
saveas(gcf,'Module_rg.png');
%% RGG
clearvars SpeactralQ;
[newman_rgg, newman_modul_rgg]= newman_comm_fast(rggA);
[newman_max_rgg,newman_max_rgg1]=max(newman_modul_rgg);
%figure;
%plot(newman_modul_reg);
%title('Newman Module');
figure;
PlotGraph(full(rggA),clustercell(newman_rgg,180,newman_max_rgg1));
title('Newman method communities RGG');
saveas(gcf,'Newman_rgg.png');
%%
rgg_spectral=GCSpectralClust2(rggA,180,8);
spectralQ=zeros(180,1);
for i=1:180
    spectralQ(i)=QFModul(rgg_spectral(:,i),rggA);
end
[spectral_max_rgg,spectral_max_rgg1]=max(spectralQ);
figure;
PlotGraph(full(rggA),rgg_spectral(:,spectral_max_rgg1));
title('Spectral method communities RGG');
saveas(gcf,'Spectral_RGG.png');
%%
Modulemax_rgg=GCModulMax1(rggA);
modulemax_val_rgg= QFModul(Modulemax_rgg,rggA);
figure;
PlotGraph(full(rggA),Modulemax_rgg);
title('Max Module method communities RGG');
saveas(gcf,'Module_rgg.png');
%% SF
clearvars SpeactralQ;
[newman_sf, newman_modul_sf]= newman_comm_fast(sf);
[newman_max_sf,newman_max_sf1]=max(newman_modul_sf);
%figure;
%plot(newman_modul_reg);
%title('Newman Module');
figure;
PlotGraph(full(sf),clustercell(newman_sf,180,newman_max_sf1));
title('Newman method communities SF');
saveas(gcf,'Newman_sf.png');
sf_spectral=GCSpectralClust2(sf,180,4);
spectralQ=zeros(180,1);
for i=1:180
    spectralQ(i)=QFModul(sf_spectral(:,i),sf);
end
[spectral_max_sf,spectral_max_sf1]=max(spectralQ);
figure;
PlotGraph(full(sf),sf_spectral(:,spectral_max_sf1));
title('Spectral method communities SF');
saveas(gcf,'Spectral_SF.png');
Modulemax_sf=GCModulMax1(sf);
modulemax_val_sf= QFModul(Modulemax_sf,sf);
figure;
PlotGraph(full(sf),Modulemax_sf);
title('Max Module method communities SF');
saveas(gcf,'Module_sf.png');
%% SW
clearvars SpeactralQ;
[newman_sw, newman_modul_sw]= newman_comm_fast(sw);
[newman_max_sw,newman_max_sw1]=max(newman_modul_sw);
%figure;
%plot(newman_modul_reg);
%title('Newman Module');
figure;
PlotGraph(full(sw),clustercell(newman_sw,180,newman_max_sw1));
title('Newman method communities SW');
saveas(gcf,'Newman_sw.png');
sw_spectral=GCSpectralClust2(sw,180,4);
spectralQ=zeros(180,1);
for i=1:180
    spectralQ(i)=QFModul(sw_spectral(:,i),sw);
end
[spectral_max_sw,spectral_max_sw1]=max(spectralQ);
figure;
PlotGraph(full(sw),sw_spectral(:,spectral_max_sw1));
title('Spectral method communities SW');
saveas(gcf,'Spectral_SW.png');
Modulemax_sw=GCModulMax1(sw);
modulemax_val_sw= QFModul(Modulemax_sw,sw);
figure;
PlotGraph(full(sw),Modulemax_sw);
title('Max Module method communities SW');
saveas(gcf,'Module_sw.png');
%% football
clearvars SpeactralQ;
[newman_foot, newman_modul_foot]= newman_comm_fast(football);
[newman_max_foot,newman_max_foot1]=max(newman_modul_foot);
%figure;
%plot(newman_modul_reg);
%title('Newman Module');
figure;
PlotGraph(full(football),clustercell(newman_foot,115,newman_max_foot1));
title('Newman method communities football');
saveas(gcf,'Newman_foot.png');
foot_spectral=GCSpectralClust2(football,115,4);
spectralQ=zeros(115,1);
for i=1:115
    spectralQ(i)=QFModul(foot_spectral(:,i),football);
end
[spectral_max_foot,spectral_max_foot1]=max(spectralQ);
figure;
PlotGraph(full(football),foot_spectral(:,spectral_max_foot1));
title('Spectral method communities football');
saveas(gcf,'Spectral_foot.png');
Modulemax_foot=GCModulMax1(football);
modulemax_val_foot= QFModul(Modulemax_foot,football);
figure;
PlotGraph(full(football),Modulemax_foot);
title('Max Module method communities football');
saveas(gcf,'Module_foot.png');
%%
%% Dolphin
clearvars SpeactralQ;
[newman_dolphins, newman_modul_dolphins]= newman_comm_fast(dolphins);
[newman_max_dolphins,newman_max_dolphins1]=max(newman_modul_dolphins);
%figure;
%plot(newman_modul_reg);
%title('Newman Module');
figure;
PlotGraph(full(dolphins),clustercell(newman_dolphins,62,newman_max_dolphins1));
title('Newman method communities dolphins');
saveas(gcf,'Newman_dolphins.png');
dolphins_spectral=GCSpectralClust2(dolphins,62,4);
spectralQ=zeros(62,1);
for i=1:62
    spectralQ(i)=QFModul(dolphins_spectral(:,i),dolphins);
end
[spectral_max_dolphins,spectral_max_dolphins1]=max(spectralQ);
figure;
PlotGraph(full(dolphins),dolphins_spectral(:,spectral_max_dolphins1));
title('Spectral method communities dolphins');
saveas(gcf,'Spectral_dolphins.png');
Modulemax_dolphins=GCModulMax1(dolphins);
modulemax_val_dolphins= QFModul(Modulemax_dolphins,dolphins);
figure;
PlotGraph(full(dolphins),Modulemax_dolphins);
title('Max Module method communities dolphins');
saveas(gcf,'Module_dolphins.png');
%% lesmis
clearvars SpeactralQ;
[newman_lesmis, newman_modul_lesmis]= newman_comm_fast(lesmis);
[newman_max_lesmis,newman_max_lesmis1]=max(newman_modul_lesmis);
%figure;
%plot(newman_modul_reg);
%title('Newman Module');
figure;
PlotGraph(full(lesmis),clustercell(newman_lesmis,77,newman_max_lesmis1));
title('Newman method communities lesmis');
saveas(gcf,'Newman_lesmis.png');
lesmis_spectral=GCSpectralClust2(lesmis,77,4);
spectralQ=zeros(77,1);
for i=1:77
    spectralQ(i)=QFModul(lesmis_spectral(:,i),lesmis);
end
[spectral_max_lesmis,spectral_max_lesmis1]=max(spectralQ);
figure;
PlotGraph(full(lesmis),lesmis_spectral(:,spectral_max_lesmis1));
title('Spectral method communities lesmis');
saveas(gcf,'Spectral_lesmis.png');
Modulemax_lesmis=GCModulMax1(lesmis);
modulemax_val_lesmis= QFModul(Modulemax_lesmis,lesmis);
figure;
PlotGraph(full(lesmis),Modulemax_lesmis);
title('Max Module method communities lesmis');
saveas(gcf,'Module_lesmis.png');
%%
[x,y]= getNodeCoordinates(115);
regNodes=[x y];
gplot(football,regNodes);