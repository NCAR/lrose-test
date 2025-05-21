% Read and diplay radar data

clear all;
close all;

addpath(genpath('~/git/lrose-test/bomb_snowstorm/analysis/'));

showPlot='on';
radar='KFTG';

figdir='/scr/cirrus1/rsfdata/projects/bomb_snowstorm/figures/vradRegPaper/';


%% Infiles

regFile='/scr/cirrus1/rsfdata/projects/nexrad/tables/DOPPLER20200525_051615VD6.txt';
lev2File='/scr/cirrus1/rsfdata/projects/nexrad/cfradial/nexrad.level2/kddc/20200525/cfrad.20200525_051303.384_to_20200525_051904.931_KDDC_SUR.nc';

%% Read Data

disp('Loading data ...')
dataReg=readDataTables(regFile,' ');

dataLev2=[];
dataLev2.VEL=[];
dataLev2.PURPLE_HAZE=[];
dataLev2=read_spol(lev2File,dataLev2);
dataLev2=dataLev2(9);

% Cut range
dataLev2.VEL(:,dataLev2.range>232)=nan;
dataLev2.PURPLE_HAZE(:,dataLev2.range>232)=nan;


nyquist=ncread(lev2File,'nyquist_velocity');
dataLev2.VEL(dataLev2.PURPLE_HAZE==1)=-99;

%% Censor regression

kernel=[9,5]; % Az and range of std kernel. Default: [9,5]
[stdVel,~]=fast_nd_std(dataReg.VEL_F,kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',mode(nyquist));

regVelCensored=dataReg.VEL_F;
regVelCensored(stdVel>9)=nan;

%% Plot preparation

xlimits1=[-235,235];
ylimits1=[-235,235];

%% Plot

close all
f1 = figure('Position',[200 500 600 1000],'DefaultAxesFontSize',12,'visible',showPlot);
colLims=[-inf,-40,-27,-21,-17,-13,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,13,17,21,27,40,inf];

t = tiledlayout(2,1,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);

ang_p = deg2rad(90-dataLev2.azimuth);

angMat=repmat(ang_p,size(dataLev2.range,1),1);

XX = (dataLev2.range.*cos(angMat));
YY = (dataLev2.range.*sin(angMat));

h1=surf(XX,YY,dataLev2.VEL,'edgecolor','none');
view(2);
title('(a) Level-2 velocity (m s^{-1})');
xlabel('km');
ylabel('km');

yticks(-300:100:300);

grid on
box on

applyColorScale(h1,dataLev2.VEL,vel_default2,colLims);

xlim(xlimits1)
ylim(ylimits1)

ar = annotation("textarrow",[0.44,0.5],[0.9,0.865],'String','1st to 2nd boundary','TextBackgroundColor','w', ...
    'TextEdgeColor','k','LineWidth',2,'TextLineWidth',1);

s2=nexttile(2);

ang_p = deg2rad(90-dataReg.azimuth);

angMat=repmat(ang_p,size(dataReg.range,1),1);

XX = (dataReg.range.*cos(angMat));
YY = (dataReg.range.*sin(angMat));

h1=surf(XX,YY,dataReg.VEL_F,'edgecolor','none');
view(2);
title('(b) REG velocity (m s^{-1})');
xlabel('km');
ylabel('km');

yticks(-300:100:300);

grid on
box on

applyColorScale(h1,dataReg.VEL_F,vel_default2,colLims);

xlim(xlimits1)
ylim(ylimits1)

ar = annotation("textarrow",[0.44,0.5],[0.4,0.36], ... 
    'String','1st to 2nd boundary','TextBackgroundColor','w','TextEdgeColor','k', ...
    'Color','k','TextColor','k','LineWidth',2,'TextLineWidth',1);

daspect(s2,[1 1 1]);
daspect(s1,[1 1 1]);

print([figdir,'figure1.png'],'-dpng','-r0');