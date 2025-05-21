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

xlimits1=[-200,200];
ylimits1=[-200,200];

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
surf(XX,YY,dataWN.DBZ,'edgecolor','none');
view(2);
clim([-3 63])
s1.Colormap=dbz_default3;
cb1=colorbar('XTick',-3:4:67);
title('(a) Reflectivity Level2 (dBZ)')
ylabel('km');

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s1,[1 1 1]);

scatter(0,0,60,'filled','MarkerFaceColor','w','MarkerEdgeColor','k');
text(-27,28,['KFTG'],'Color','k','FontSize',12,'FontWeight','bold');

s1.SortMethod='childorder';

% CMD

s2=nexttile(2);
h=surf(XX2,YY2,dataR.CMD_FLAG,'edgecolor','none');
view(2);
title('(b) CMD flag')

s2.Colormap=[0,0,1;1,0,0];
clim([0,1]);
colorbar('Ticks',[0.25,0.75],'TickLabels',{'0','1'});

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s2,[1 1 1]);


% Velocity

s3=nexttile(3);
h2=surf(XX,YY,dataWN.VEL,'edgecolor','none');
view(2);
title('(c) Velocity Level2 WN (m s^{-1})')
xlabel('km');
ylabel('km');

grid on
box on

colLims=[-inf,-30,-26,-21,-17,-13,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,13,17,21,26,30,inf];
applyColorScale(h2,dataWN.VEL,vel_default2,colLims);

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s3,[1 1 1]);

% ORDER

s4=nexttile(4);
dataR.REGR_ORDER(dataR.CMD_FLAG==0)=nan;
h=surf(XX2,YY2,dataR.REGR_ORDER,'edgecolor','none');
view(2);
title('(d) Polynomial order')
xlabel('km');

orderMax=max(dataR.REGR_ORDER(:),[],'omitmissing');
orderMin=min(dataR.REGR_ORDER(:),[],'omitmissing');
s4.Colormap=cat(1,turbo(orderMax-orderMin+1));
clim([orderMin-0.5,orderMax+0.5]);
colorbar('Ticks',2:11,'TickLabels',{'2','3','4','5','6','7','8','9','10','11'});

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s4,[1 1 1]);

%print([figdir,'figure21.png'],'-dpng','-r0');
exportgraphics(f1,[figdir,'figure18.png'],'Resolution','300');