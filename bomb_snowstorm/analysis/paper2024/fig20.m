% Read and diplay radar data

clear all;
close all;

addpath(genpath('~/git/lrose-test/bomb_snowstorm/analysis/'));

showPlot='on';

figdir='/scr/cirrus1/rsfdata/projects/bomb_snowstorm/figures/paper2024/';

%% Load data

infileWN='/scr/sleet1/rsfdata/projects/eolbase/cfradial/kftg/moments/20220329/cfrad.20220329_221646.829_to_20220329_222242.984_KFTG_SUR.nc';

dataWN=[];

dataWN.DBZ=[];
dataWN.VEL=[];
dataWN.WIDTH=[];
dataWN.ZDR=[];
dataWN.PHIDP=[];
dataWN.RHOHV=[];

dataWN=read_spol(infileWN,dataWN);

dataWN=dataWN(1);

infileR='/scr/cirrus1/rsfdata/projects/nexrad/tables/KFTG_LPRT_SR_20220329_221645_1.91_253.46_14ptFlt28pt-V4.txt';
dataR=readDataTables(infileR,' ');

%% Plot preparation
ang_p = deg2rad(90-dataWN.azimuth);

angMat=repmat(ang_p,size(dataWN.range,1),1);

xlimits1=[-120,60];
ylimits1=[-110,70];

XX = (dataWN.range.*cos(angMat));
YY = (dataWN.range.*sin(angMat));

ang_p2 = deg2rad(90-dataR.azimuth);

angMat2=repmat(ang_p2,size(dataR.range,1),1);

XX2 = (dataR.range.*cos(angMat2));
YY2 = (dataR.range.*sin(angMat2));

%% Plot
close all

f1=figure('Position',[200 500 1000 900],'DefaultAxesFontSize',12,'visible',showPlot);
t = tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);
hold on
h1=surf(XX,YY,dataWN.PHIDP,'edgecolor','none');
view(2);
clim([0,114]);
title(['(a) \phi_{DP} Level2 (',char(176),')'])
ylabel('km');
s1.Colormap=phidp_default;
cb1=colorbar;
cb1.Ticks=0:9:120;

grid on
box on

scatter(0,0,60,'filled','MarkerFaceColor','w','MarkerEdgeColor','k');
text(-10,7,['KFTG'],'Color','w','FontSize',12,'FontWeight','bold');

xlim(xlimits1)
ylim(ylimits1)
daspect(s1,[1 1 1]);

% rectangle('Position',[5 -17 40 55],'EdgeColor','w','LineWidth',1.5);
% scatter(0,0,90,'filled','MarkerFaceColor','w','MarkerEdgeColor','k');
% text(-20,0,['S-Pol'],'Color','w','FontSize',12,'FontWeight','bold');

s1.SortMethod='childorder';

% Width Reg.

s2=nexttile(2);
hold on
h3=surf(XX2,YY2,dataR.PHIDP_F,'edgecolor','none');
view(2);
clim([0,114]);
title(['(b) \phi_{DP} Regression (',char(176),')'])
s2.Colormap=phidp_default;
cb3=colorbar;
cb3.Ticks=0:9:120;

grid on
box on

scatter(0,0,60,'filled','MarkerFaceColor','w','MarkerEdgeColor','k');

xlim(xlimits1)
ylim(ylimits1)
daspect(s2,[1 1 1]);

%rectangle('Position',[5 -17 40 55],'EdgeColor','w','LineWidth',1.5);

s2.SortMethod='childorder';

% Rho WN

s3=nexttile(3);
hold on
h2=surf(XX,YY,dataWN.RHOHV,'edgecolor','none');
view(2);
title('(c) \rho_{HV} Level2')
xlabel('km');
ylabel('km');

scatter(0,0,60,'filled','MarkerFaceColor','w','MarkerEdgeColor','k');

grid on
box on

colLims=[-inf,0,0.7,0.8,0.85,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.975,0.98,0.985,0.99,0.995,1.1,inf];
applyColorScale(h2,dataWN.RHOHV,rhohv_default,colLims);

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s3,[1 1 1]);

%rectangle('Position',[5 -17 40 55],'EdgeColor','w','LineWidth',1.5);

s3.SortMethod='childorder';

% ZDR Reg.

s4=nexttile(4);
hold on
h4=surf(XX2,YY2,dataR.RHOHV_NNC_F,'edgecolor','none');
view(2);
title('(d) \rho_{HV} Regression')
xlabel('km');

grid on
box on

colLims=[-inf,0,0.7,0.8,0.85,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.975,0.98,0.985,0.99,0.995,1.1,inf];
applyColorScale(h4,dataR.RHOHV_NNC_F,rhohv_default,colLims);

grid on
box on

scatter(0,0,60,'filled','MarkerFaceColor','w','MarkerEdgeColor','k');

xlim(xlimits1)
ylim(ylimits1)
daspect(s4,[1 1 1]);

%rectangle('Position',[5 -17 40 55],'EdgeColor','w','LineWidth',1.5);

s4.SortMethod='childorder';

%print([figdir,'figure23.png'],'-dpng','-r0');
exportgraphics(f1,[figdir,'figure20.png'],'Resolution','300');