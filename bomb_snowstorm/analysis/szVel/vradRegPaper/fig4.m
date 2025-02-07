% Read and diplay radar data

clear all;
close all;

addpath(genpath('~/git/lrose-test/bomb_snowstorm/analysis/'));

showPlot='on';

figdir='/scr/cirrus1/rsfdata/projects/bomb_snowstorm/figures/vradRegPaper/';

%% Load data

infileWN='/scr/cirrus1/rsfdata/projects/nexrad/cfradial/nexrad.level2/kftg/20220329/cfrad.20220329_190500.493_to_20220329_191029.825_KFTG_SUR.nc';
dataWN.REF=[];
dataWN=read_spol(infileWN,dataWN);

dataWN=dataWN(1);

%% Plot preparation

ang_p = deg2rad(90-dataWN.azimuth);

angMat=repmat(ang_p,size(dataWN.range,1),1);
xlimits1=[-200,260];
ylimits1=[-220,220];

XX = (dataWN.range.*cos(angMat));
YY = (dataWN.range.*sin(angMat));

%% Plot
close all

figure('Position',[200 500 600 510],'DefaultAxesFontSize',12,'visible',showPlot);
t = tiledlayout(1,1,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);
hold on
surf(XX,YY,dataWN.REF,'edgecolor','none');
view(2);
clim([-3 63])
s1.Colormap=dbz_default3;
cb1=colorbar('XTick',-3:4:67);
xlabel('km');
ylabel('km');

xlim(xlimits1);
ylim(ylimits1);

daspect(s1,[1 1 1]);

title('Reflectivity (dBZ)')

xlabel('km');
ylabel('km');

grid on
box on

%scatter(0,0,90,'filled','MarkerFaceColor','w','MarkerEdgeColor','k');
%text(-15,0,['S-Pol'],'Color','w','FontSize',12,'FontWeight','bold');

%rectangle('Position',[5 -17 40 55],'EdgeColor','w','LineWidth',1.5);

s1.SortMethod='childorder';

print([figdir,'figure4.png'],'-dpng','-r0');
