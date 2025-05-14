% Read and diplay radar data

clear all;
close all;

addpath(genpath('~/git/lrose-test/bomb_snowstorm/analysis/'));

showPlot='on';

figdir='/scr/cirrus1/rsfdata/projects/bomb_snowstorm/figures/paper2024/';

%% Load data

infileWN='/scr/cirrus1/rsfdata/projects/bomb_snowstorm/tables/SPOL20190313_220622_INDX_CMD_RHV_GAUSS_WN_V3_1.txt';
dataWN=readDataTables(infileWN,' ');

infileR='/scr/cirrus1/rsfdata/projects/bomb_snowstorm/tables/SPOL20190313_220622_INDX_CMD_RHV_GAUSS_REG_V3.txt';
dataR=readDataTables(infileR,' ');

%% Plot preparation

ang_p = deg2rad(90-dataWN.azimuth);

angMat=repmat(ang_p,size(dataWN.range,1),1);

xlimits1=[-20,100];
ylimits1=[-20,100];

XX = (dataWN.range.*cos(angMat));
YY = (dataWN.range.*sin(angMat));

%% Plot
close all

f1=figure('Position',[200 500 1000 900],'DefaultAxesFontSize',12,'visible',showPlot);
t = tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);
hold on
h1=surf(XX,YY,dataWN.WIDTH_F,'edgecolor','none');
view(2);
clim([-3 47])
title('(a) Spectrum width WN (m s^{-1})')
ylabel('km');

colLims=[-inf,0,0.5,1,1.5,2,2.5,3,4,5,6,7,8,10,12.5,15,20,25,50,inf];
applyColorScale(h1,dataWN.WIDTH_F,width_default,colLims);

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s1,[1 1 1]);

rectangle('Position',[5 -17 40 55],'EdgeColor','w','LineWidth',1.5);
scatter(0,0,90,'filled','MarkerFaceColor','w','MarkerEdgeColor','k');
text(-20,2,['S-Pol'],'Color','k','FontSize',12,'FontWeight','bold');

s1.SortMethod='childorder';

% Width Reg.

s2=nexttile(2);
hold on
h3=surf(XX,YY,dataR.WIDTH_F,'edgecolor','none');
view(2);
clim([-3 47])
title('(b) Spectrum width Regression (m s^{-1})')

colLims=[-inf,0,0.5,1,1.5,2,2.5,3,4,5,6,7,8,10,12.5,15,20,25,50,inf];
applyColorScale(h3,dataR.WIDTH_F,width_default,colLims);

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s2,[1 1 1]);

rectangle('Position',[5 -17 40 55],'EdgeColor','w','LineWidth',1.5);

s2.SortMethod='childorder';


% ZDR WN

s3=nexttile(3);
h2=surf(XX,YY,dataWN.ZDR_F,'edgecolor','none');
view(2);
title('(c) Z_{DR} WN (dB)')
xlabel('km');
ylabel('km');

grid on
box on

colLims=[-inf,-20,-2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.5,2,2.5,3,4,5,6,8,10,15,20,50,99,inf];
applyColorScale(h2,dataWN.ZDR_F,zdr_default,colLims);

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s3,[1 1 1]);

rectangle('Position',[5 -17 40 55],'EdgeColor','w','LineWidth',1.5);

s3.SortMethod='childorder';

% ZDR Reg.

s4=nexttile(4);
h4=surf(XX,YY,dataR.ZDR_F,'edgecolor','none');
view(2);
title('(d) Z_{DR} Regression (dB)')
xlabel('km');

grid on
box on

colLims=[-inf,-20,-2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.5,2,2.5,3,4,5,6,8,10,15,20,50,99,inf];
applyColorScale(h4,dataR.ZDR_F,zdr_default,colLims);

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s4,[1 1 1]);

rectangle('Position',[5 -17 40 55],'EdgeColor','w','LineWidth',1.5);

s4.SortMethod='childorder';

%print([figdir,'figure14.png'],'-dpng','-r0');
exportgraphics(f1,[figdir,'figure11.png'],'Resolution','300');