% Read and diplay radar data

clear all;
close all;

addpath(genpath('~/git/lrose-test/bomb_snowstorm/analysis/'));

showPlot='on';

figdir='/scr/cirrus1/rsfdata/projects/bomb_snowstorm/figures/paper2024/';

%% Load data

infileWN='/scr/sleet1/rsfdata/projects/eolbase/cfradial/kftg/moments/20220329/cfrad.20220329_221646.829_to_20220329_222242.984_KFTG_SUR.nc';
dataWN=[];
dataWN.DBZ_F=[];
dataWN.VEL_F=[];
dataWN.WIDTH_F=[];
dataWN.ZDR_F=[];
dataWN.PHIDP_F=[];
dataWN.RHOHV_NNC_F=[];
dataWN.REGR_ORDER=[];
dataWN.CMD_FLAG=[];

dataWN=read_spol(infileWN,dataWN);

infileR='/scr/cirrus1/rsfdata/projects/bomb_snowstorm/tables/SPOL20190313_220622_INDX_CMD_RHV_GAUSS_REG_V3.txt';
dataR=readDataTables(infileR,' ');

%% Plot preparation

ang_p = deg2rad(90-dataWN.azimuth);

angMat=repmat(ang_p,size(dataWN.range,1),1);

xlimits1=[-250,250];
ylimits1=[-250,250];

XX = (dataWN.range.*cos(angMat));
YY = (dataWN.range.*sin(angMat));

%% Plot
close all

figure('Position',[200 500 1000 900],'DefaultAxesFontSize',12,'visible',showPlot);
t = tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);
hold on
surf(XX,YY,dataUF.DBZ_F,'edgecolor','none');
view(2);
clim([-10 65])
title('(a) Reflectivity (dBZ)')
ylabel('km');
s1.Colormap=dbz_default2;
cb1=colorbar('XTick',-10:3:65);

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

% Notch width

s2=nexttile(2);
dataWN.REGR_ORDER(dataR.CMD_FLAG==0)=nan;
dataWN.REGR_ORDER(dataWN.REGR_ORDER==0)=nan;
h=surf(XX,YY,dataWN.REGR_ORDER,'edgecolor','none');
view(2);
title('(c) Notch width')

orderMax=max(dataWN.REGR_ORDER(:),[],'omitmissing');
orderMin=min(dataWN.REGR_ORDER(:),[],'omitmissing');
s2.Colormap=turbo(orderMax-orderMin+1);
clim([orderMin-0.5,orderMax+0.5]);
colorbar('Ticks',5:17,'TickLabels',{'5','6','7','8','9','10','11','12','13','14','15','16','17'});

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s2,[1 1 1]);

% CMD

s3=nexttile(3);
h=surf(XX,YY,dataR.CMD_FLAG,'edgecolor','none');
view(2);
title('(b) CMD flag')
xlabel('km');
ylabel('km');

s3.Colormap=[0,0,1;1,0,0];
clim([0,1]);
colorbar('Ticks',[0.25,0.75],'TickLabels',{'0','1'});

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

print([figdir,'figure22.png'],'-dpng','-r0');
