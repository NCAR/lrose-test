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
surf(XX,YY,dataWN.DBZ_F,'edgecolor','none');
view(2);
clim([-3 63])
s1.Colormap=dbz_default3;
cb1=colorbar('XTick',-3:4:67);
title('(a) Reflectivity WN (dBZ)')
ylabel('km');

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s1,[1 1 1]);

rectangle('Position',[5 -17 40 55],'EdgeColor','w','LineWidth',1.5);
scatter(0,0,90,'filled','MarkerFaceColor','w','MarkerEdgeColor','k');
text(-20,0,['S-Pol'],'Color','w','FontSize',12,'FontWeight','bold');

s1.SortMethod='childorder';

% Refl. Reg.

s2=nexttile(2);
hold on
surf(XX,YY,dataR.DBZ_F,'edgecolor','none');
view(2);
clim([-3 63])
s2.Colormap=dbz_default3;
cb1=colorbar('XTick',-3:4:67);
title('(b) Reflectivity Regression (dBZ)')

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s2,[1 1 1]);

rectangle('Position',[5 -17 40 55],'EdgeColor','w','LineWidth',1.5);

s2.SortMethod='childorder';

% Vel. WN

s3=nexttile(3);
h2=surf(XX,YY,dataWN.VEL_F,'edgecolor','none');
view(2);
title('(c) Velocity WN (m s^{-1})')
xlabel('km');
ylabel('km');


grid on
box on

colLims=[-inf,-30,-26,-21,-17,-13,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,13,17,21,26,30,inf];
applyColorScale(h2,dataWN.VEL_F,vel_default2,colLims);

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s3,[1 1 1]);

rectangle('Position',[5 -17 40 55],'EdgeColor','w','LineWidth',1.5);

s3.SortMethod='childorder';

% Vel Reg.

s4=nexttile(4);
h4=surf(XX,YY,dataR.VEL_F,'edgecolor','none');
view(2);
title('(d) Velocity Regression (m s^{-1})')
xlabel('km');

grid on
box on

colLims=[-inf,-30,-26,-21,-17,-13,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,13,17,21,26,30,inf];
applyColorScale(h4,dataWN.VEL_F,vel_default2,colLims);

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s4,[1 1 1]);

rectangle('Position',[5 -17 40 55],'EdgeColor','w','LineWidth',1.5);

s4.SortMethod='childorder';

%print([figdir,'figure13.png'],'-dpng','-r0');
exportgraphics(f1,[figdir,'figure10.png'],'Resolution','300');