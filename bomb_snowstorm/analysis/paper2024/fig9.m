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

infileUF='/scr/cirrus1/rsfdata/projects/bomb_snowstorm/tables/SPOL20190313_220622_INDX_UNFILTERED.txt';
dataUF=readDataTables(infileUF,' ');

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
surf(XX,YY,dataUF.DBZ_F,'edgecolor','none');
view(2);
clim([-3 63])
s1.Colormap=dbz_default3;
cb1=colorbar('XTick',-3:4:67);
title('(a) Reflectivity (dBZ)')
ylabel('km');

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s1,[1 1 1]);

rectangle('Position',[-18 -18 17 98],'EdgeColor','w','LineWidth',1.5);
text(2,73,[{'Rocky'};{'Mountains'}],'Color','w','FontSize',12,'FontWeight','bold');
scatter(0,0,90,'filled','MarkerFaceColor','w','MarkerEdgeColor','k');
text(5,0,['S-Pol'],'Color','w','FontSize',16,'FontWeight','bold');

s1.SortMethod='childorder';

% CMD

s2=nexttile(2);
h=surf(XX,YY,dataR.CMD_FLAG,'edgecolor','none');
view(2);
title('(b) CMD flag')

s2.Colormap=[0,0,1;1,0,0];
clim([0,1]);
colorbar('Ticks',[0.25,0.75],'TickLabels',{'0','1'});

rectangle('Position',[5 -17 40 55],'EdgeColor','w','LineWidth',1.5);
s2.SortMethod='childorder';

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s2,[1 1 1]);


% Notch width

s3=nexttile(3);
dataWN.REGR_ORDER(dataR.CMD_FLAG==0)=nan;
dataWN.REGR_ORDER(dataWN.REGR_ORDER==0)=nan;
h=surf(XX,YY,dataWN.REGR_ORDER,'edgecolor','none');
view(2);
title('(c) Notch width')
xlabel('km');
ylabel('km');

orderMax=max(dataWN.REGR_ORDER(:),[],'omitmissing');
orderMin=min(dataWN.REGR_ORDER(:),[],'omitmissing');
s3.Colormap=turbo(orderMax-orderMin+1);
clim([orderMin-0.5,orderMax+0.5]);
colorbar('Ticks',5:17,'TickLabels',{'5','6','7','8','9','10','11','12','13','14','15','16','17'});

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s3,[1 1 1]);

% ORDER

s4=nexttile(4);
dataR.REGR_ORDER(dataR.CMD_FLAG==0)=nan;
h=surf(XX,YY,dataR.REGR_ORDER,'edgecolor','none');
view(2);
title('(d) Polynomial order')
xlabel('km');

orderMax=12;
orderMin=min(dataR.REGR_ORDER(:),[],'omitmissing');
s4.Colormap=cat(1,turbo(orderMax-orderMin+1),[1,0,1]);
clim([orderMin-0.5,orderMax+1.5]);
colorbar('Ticks',3:13,'TickLabels',{'3','4','5','6','7','8','9','10','11','12','21'});

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s4,[1 1 1]);

%print([figdir,'figure12.png'],'-dpng','-r0');
exportgraphics(f1,[figdir,'figure9.png'],'Resolution','300');