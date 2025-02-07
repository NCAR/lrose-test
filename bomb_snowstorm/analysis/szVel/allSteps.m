% Read and diplay radar data

clear all;
close all;

addpath(genpath('~/git/lrose-test/bomb_snowstorm/analysis/'));

showPlot='on';
radar='KFTG';
%radar='KLWX';
%radar='KDDC1';
%radar='KDDC2';

figdir=['/scr/cirrus1/rsfdata/projects/nexrad/figures/szComp/',radar,'/'];

nyquist=25.7428;

%% Infiles
if strcmp(radar,'KFTG')
    regFile='/scr/cirrus1/rsfdata/projects/nexrad/tables/KFTG_SZ_20220329_190532_0.48_272.92_O37_64pts_V5.txt';
    regVradFile='/scr/cirrus1/rsfdata/projects/nexrad/matFiles/KFTG_Regression_and_VRAD_Filt_13.mat';
    xlimits1=[-200,260];
    ylimits1=[-220,220];
    xlimits2=[50,150];
    ylimits2=[70,170];
    censThreshStd=7;
elseif strcmp(radar,'KLWX')
    regFile='/scr/cirrus1/rsfdata/projects/nexrad/tables/DOPklwx20230807_220627.txt';
    regVradFile='/scr/cirrus1/rsfdata/projects/nexrad/matFiles/KLWX_Regression_and_VRAD_Filt_12.mat';
    xlimits1=[-260,200];
    ylimits1=[-170,270];
    xlimits2=[50,170];
    ylimits2=[60,180];
    censThreshStd=7;
elseif strcmp(radar,'KDDC1')
    regFile='/scr/cirrus1/rsfdata/projects/nexrad/tables/DOPPLER20200525_051615VD6.txt';
    regVradFile='/scr/cirrus1/rsfdata/projects/nexrad/matFiles/KDDC_Scan1_Regression_and_VRAD.mat';
    xlimits1=[-230,230];
    ylimits1=[-230,230];
    xlimits2=[-150,150];
    ylimits2=[-150,150];
    censThreshStd=7;
elseif strcmp(radar,'KDDC2')
    regFile='/scr/cirrus1/rsfdata/projects/nexrad/tables/DOPPLER20200525_053516VD6.txt';
    regVradFile='/scr/cirrus1/rsfdata/projects/nexrad/matFiles/KDDC_Scan2_Regression_and_VRAD.mat';
    xlimits1=[-300,300];
    ylimits1=[-300,300];
    xlimits2=[-200,200];
    ylimits2=[-200,200];
    censThreshStd=7;
end

%% Read regressino data

disp('Loading data ...')
regIn=readDataTables(regFile,' ');

%% Read mat files

regVradIn=load(regVradFile);

%% Match azimuths

azRes=round((regIn.azimuth(2)-regIn.azimuth(1))*10)/10;
% if azRes==0.5
%     pastDot=data1in.azimuth(1)-floor(data1in.azimuth(1));
%     if (pastDot>=0.2 & pastDot<=0.3) | (pastDot>=0.7 & pastDot<=0.8)
%         allAz=0.25:azRes:360;
%     else
%         allAz=0.5:azRes:360;
%     end
% else
    allAz=0:360;
% end

ib1=[];
ib2=[];
ibAll=[];
for kk=1:length(allAz)
    [minDiff1,minInd1]=min(abs(regIn.azimuth-allAz(kk)));
    [minDiff2,minInd2]=min(abs(regVradIn.az-allAz(kk)));
    %if minDiff1<azRes/2+0.01 & minDiff2<azRes/2+0.01
        ib1=cat(1,ib1,minInd1);
        ib2=cat(1,ib2,minInd2);
        ibAll=cat(1,ibAll,kk);
    %end
end

reg=[];
reg.range=regIn.range;
regVrad=[];
lev2=[];
vradLeg=[];

% Regression
reg.VEL=nan(length(allAz),size(regIn.VEL_F,2));
reg.VEL(ibAll,:)=regIn.VEL_F(ib1,:);
% Regression and vrad
vradT=(regVradIn.vrad.*regVradIn.thr)';
regVrad.VEL=nan(length(allAz),size(vradT,2));
regVrad.VEL(ibAll,:)=vradT(ib2,:);
% Regression input for VRAD
vradInT=(regVradIn.v1)';
regInVrad.VEL=nan(length(allAz),size(vradInT,2));
regInVrad.VEL(ibAll,:)=vradInT(ib2,:);

%% Cut off range
maxRange=min([regIn.range(end),regVradIn.range_km(end)]);
minRange=max([regIn.range(1),regVradIn.range_km(1)]);

ri1=find(reg.range>=minRange & reg.range<=maxRange);
reg.VEL=reg.VEL(:,ri1);
reg.range=reg.range(:,ri1);

ri2=find(regVradIn.range_km>=minRange & regVradIn.range_km<=maxRange);
regVrad.VEL=regVrad.VEL(:,ri2);
regInVrad.VEL=regInVrad.VEL(:,ri2);


%% Censor regression

regCensored=reg.VEL;

kernel=[9,5]; % Az and range of std kernel. Default: [9,5]
[stdVel,~]=fast_nd_std(reg.VEL,kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',mode(nyquist));
regCensored(stdVel>censThreshStd)=nan;

%% Fill in with regression
regVradFilled=regVrad.VEL;
regVradFilled(isnan(regVradFilled))=regCensored(isnan(regVradFilled));

%% Plot preparation

ang_p = deg2rad(90-allAz);

angMat=repmat(ang_p',size(reg.range,1),1);

XX = (reg.range.*cos(angMat));
YY = (reg.range.*sin(angMat));

%% Plot

close all
f1 = figure('Position',[200 500 1200 1000],'DefaultAxesFontSize',12,'visible',showPlot);
colLims=[-inf,-40,-27,-21,-17,-13,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,13,17,21,27,40,inf];

t = tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

% Regression original
s1=nexttile(1);

h1=surf(XX,YY,reg.VEL,'edgecolor','none');
view(2);
title('VEL regression (m s^{-1})');
xlabel('km');
ylabel('km');

grid on
box on

applyColorScale(h1,reg.VEL,vel_default2,colLims);

xlim(xlimits1)
ylim(ylimits1)
daspect(s1,[1 1 1]);

% Regression censored
s2=nexttile(2);
regCensoredVRAD=regInVrad.VEL;
h1=surf(XX,YY,regCensoredVRAD,'edgecolor','none');
view(2);
title('VEL regression VRAD input (m s^{-1})');
xlabel('km');
ylabel('km');

grid on
box on

applyColorScale(h1,regCensoredVRAD,vel_default2,colLims);

xlim(xlimits1)
ylim(ylimits1)
daspect(s2,[1 1 1]);

% VRAD and regression
s3=nexttile(3);

h1=surf(XX,YY,regVrad.VEL,'edgecolor','none');
view(2);
title('VEL VRAD and regression (m s^{-1})');
xlabel('km');
ylabel('km');

grid on
box on

applyColorScale(h1,regVrad.VEL,vel_default2,colLims);

xlim(xlimits1)
ylim(ylimits1)
daspect(s3,[1 1 1]);

% VRAD and regression filled
s4=nexttile(4);

h1=surf(XX,YY,regVradFilled,'edgecolor','none');
view(2);
title('VEL VRAD and regression filled (m s^{-1})');
xlabel('km');
ylabel('km');

grid on
box on

applyColorScale(h1,regVradFilled,vel_default2,colLims);

xlim(xlimits1)
ylim(ylimits1)
daspect(s4,[1 1 1]);

linkaxes([s1,s2,s3,s4],'xy');

print([figdir,radar,'_PPIs.png'],'-dpng','-r0');

% Second zoom
s1.XLim=xlimits2;
s1.YLim=ylimits2;
daspect(s1,[1 1 1]);
daspect(s2,[1 1 1]);
daspect(s3,[1 1 1]);
daspect(s4,[1 1 1]);

print([figdir,radar,'_PPIs_zoom.png'],'-dpng','-r0');

%% Area coverage

pixNum=[sum(~isnan(lev2.VEL(:))),sum(~isnan(vradLeg.VEL(:))),sum(~isnan(regVradFilled(:)))];
pixPerc=pixNum./pixNum(1).*100;

%close all
f1 = figure('Position',[200 500 600 500],'DefaultAxesFontSize',12,'visible',showPlot);

t = tiledlayout(1,1,'TileSpacing','tight','Padding','compact');
s1=nexttile(1);

bar(pixPerc,1);
xlim([0.5,3.5]);

ylabel('Percent (%)')
xticklabels({'Level 2';'Level 2 + VRAD';'Regression + VRAD'});
xtickangle(0);

title('Increase in coverage')

print([figdir,radar,'_coverageIncrease.png'],'-dpng','-r0');

%% Smoothness

vradCensLev2=vradLeg.VEL;
vradCensLev2(isnan(lev2.VEL))=nan;
vradRegCensLev2=regVradFilled;
vradRegCensLev2(isnan(lev2.VEL))=nan;
vradRegCensVrad=regVradFilled;
vradRegCensVrad(isnan(vradLeg.VEL))=nan;

[stdLev2,~]=fast_nd_std(lev2.VEL,kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',nyquist);
[stdVradLeg,~]=fast_nd_std(vradLeg.VEL,kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',nyquist);
[stdVradLegClev2,~]=fast_nd_std(vradCensLev2,kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',nyquist);
[stdVradRegVrad,~]=fast_nd_std(vradRegCensVrad,kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',nyquist);
[stdVradRegClev2,~]=fast_nd_std(vradRegCensLev2,kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',nyquist);

stdLev2(isnan(lev2.VEL))=nan;
stdVradLegClev2(isnan(lev2.VEL))=nan;
stdVradRegClev2(isnan(lev2.VEL))=nan;
stdVradLeg(isnan(vradLeg.VEL))=nan;
stdVradRegVrad(isnan(vradLeg.VEL))=nan;

% Plot
stdLims=[0,7];
diffLims=[-3,3];

% close all
f1 = figure('Position',[200 500 1500 1200],'DefaultAxesFontSize',12,'visible',showPlot);

t = tiledlayout(3,3,'TileSpacing','tight','Padding','compact');

% NEXRAD level 2
s1=nexttile(1);

h1=surf(XX,YY,stdLev2,'edgecolor','none');
view(2);
s1.Colormap=jet;
clim(stdLims);
colorbar
title('St. dev. NEXRAD level 2 (m s^{-1})');
xlabel('km');
ylabel('km');

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s1,[1 1 1]);

s2=nexttile(2);
vradMinLev2=stdVradLegClev2-stdLev2;

h1=surf(XX,YY,vradMinLev2,'edgecolor','none');
view(2);
s2.Colormap=velCols;
clim(diffLims);
colorbar
title('St. dev. VRAD/lev. 2 - St. dev. level 2 (m s^{-1})');
xlabel('km');
ylabel('km');

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s2,[1 1 1]);

s3=nexttile(3);

hold on
edges=-1.05:0.1:1.05;
hc=histcounts(vradMinLev2(:),edges);
bar(edges(1:end-1)+(edges(2)-edges(1))/2,hc/sum(hc)*100,1)
xlim([-0.5,0.5]);

ylims=s3.YLim;
plot([0,0],ylims,'-r','LineWidth',2);

xlabel('Velocity (m s^{-1})')
ylabel('Percent (%)')

s3.SortMethod='childorder';

grid on
box on
title('St. dev. VRAD/lev. 2 - St. dev. level 2 (m s^{-1})');

% VRAD
s4=nexttile(4);

h1=surf(XX,YY,stdVradLeg,'edgecolor','none');
view(2);
s4.Colormap=jet;
clim(stdLims);
colorbar
title('St. dev. VRAD/lev. 2 (m s^{-1})');
xlabel('km');
ylabel('km');

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s4,[1 1 1]);

s5=nexttile(5);
vradRegMinVrad=stdVradRegVrad-stdVradLeg;

h1=surf(XX,YY,vradRegMinVrad,'edgecolor','none');
view(2);
s5.Colormap=velCols;
clim(diffLims);
colorbar
title('St. dev. VRAD/REGR - St. dev. VRAD/lev. 2 (m s^{-1})');
xlabel('km');
ylabel('km');

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s5,[1 1 1]);

s6=nexttile(6);

hold on
edges=-4:0.08:4;
hc=histcounts(vradRegMinVrad(:),edges);
bar(edges(1:end-1)+(edges(2)-edges(1))/2,hc/sum(hc)*100,1)
xlim([-0.7,0.7]);

ylims=s6.YLim;
plot([0,0],ylims,'-r','LineWidth',2);

xlabel('Velocity (m s^{-1})')
ylabel('Percent (%)')

s6.SortMethod='childorder';

grid on
box on
title('St. dev. VRAD/REGR - St. dev. VRAD/lev. 2 (m s^{-1})');

% VRAD Regr
s7=nexttile(7);

h1=surf(XX,YY,stdVradRegVrad,'edgecolor','none');
view(2);
s7.Colormap=jet;
clim(stdLims);
colorbar
title('St. dev. VRAD/REGR (m s^{-1})');
xlabel('km');
ylabel('km');

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s7,[1 1 1]);

s8=nexttile(8);
vradRegMinLev2=stdVradRegClev2-stdLev2;

h1=surf(XX,YY,vradRegMinLev2,'edgecolor','none');
view(2);
s8.Colormap=velCols;
clim(diffLims);
colorbar
title('St. dev. VRAD/REGR - St. dev. level 2 (m s^{-1})');
xlabel('km');
ylabel('km');

grid on
box on

xlim(xlimits1)
ylim(ylimits1)
daspect(s8,[1 1 1]);

s9=nexttile(9);

hold on
edges=-4:0.08:4;
hc=histcounts(vradRegMinLev2(:),edges);
bar(edges(1:end-1)+(edges(2)-edges(1))/2,hc/sum(hc)*100,1)
xlim([-0.7,0.7]);

ylims=s9.YLim;
plot([0,0],ylims,'-r','LineWidth',2);

xlabel('Velocity (m s^{-1})')
ylabel('Percent (%)')

s9.SortMethod='childorder';

grid on
box on
title('St. dev. VRAD/REGR - St. dev. level 2 (m s^{-1})');

linkaxes([s1,s2,s4,s5,s7,s8],'xy');

print([figdir,radar,'_stDev.png'],'-dpng','-r0');

% Second zoom
s1.XLim=xlimits2;
s1.YLim=ylimits2;
daspect(s1,[1 1 1]);
daspect(s2,[1 1 1]);
daspect(s4,[1 1 1]);
daspect(s5,[1 1 1]);
daspect(s7,[1 1 1]);
daspect(s8,[1 1 1]);

print([figdir,radar,'_stDev_zoom.png'],'-dpng','-r0');

%% Histograms

phPix=isnan(lev2.VEL) & ~isnan(vradLeg.VEL);

vradPH=vradLeg.VEL(phPix==1);
regPH=reg.VEL(phPix==1);
regVradPH=regVradFilled(phPix==1);

% Plot

f1 = figure('Position',[200 500 700 1100],'DefaultAxesFontSize',12,'visible',showPlot);
t = tiledlayout(3,1,'TileSpacing','tight','Padding','compact');
s1=nexttile(1);

hold on
edges=-20.25:0.5:20.25;
hc=histcounts(regPH-vradPH,edges);
bar(edges(1:end-1)+(edges(2)-edges(1))/2,hc/sum(hc)*100,1)
xlim([-20,20]);

title('Regression minus VRAD/level 2')

xlabel('Velocity (m s^{-1})')
ylabel('Percent (%)')

grid on
box on

s2=nexttile(2);

hold on
edges=-20.25:0.5:20.25;
hc=histcounts(regVradPH-vradPH,edges);
bar(edges(1:end-1)+(edges(2)-edges(1))/2,hc/sum(hc)*100,1)
xlim([-20,20]);

title('Regression/VRAD minus VRAD/level 2')

xlabel('Velocity (m s^{-1})')
ylabel('Percent (%)')

grid on
box on

s3=nexttile(3);
hold on
edges=-20.25:0.5:20.25;
hc1=histcounts(regPH-vradPH,edges);
stairs(edges(1:end-1)+(edges(2)-edges(1))/2,hc1,'LineWidth',2)
hc2=histcounts(regVradPH-vradPH,edges);
stairs(edges(1:end-1)+(edges(2)-edges(1))/2,hc2,'LineWidth',2)
xlim([-20,20]);

legend({'Regression - VRAD/level 2','Regression/VRAD - VRAD/level 2'},'Orientation','horizontal','Location','northoutside');

xlabel('Velocity (m s^{-1})')
ylabel('Count')

grid on
box on

print([figdir,radar,'_histPurpleHaze.png'],'-dpng','-r0');
<<<<<<< HEAD
<<<<<<< HEAD

%% Plot purple haze

vradPHsurf=vradLeg.VEL;
vradPHsurf(phPix==0)=nan;
regPHsurf=reg.VEL;
regPHsurf(phPix==0)=nan;
regVradPHsurf=regVradFilled;
regVradPHsurf(phPix==0)=nan;

f1 = figure('Position',[200 500 500 1200],'DefaultAxesFontSize',12,'visible',showPlot);
t = tiledlayout(3,1,'TileSpacing','tight','Padding','compact');
s1=nexttile(1);

h1=surf(XX,YY,vradPHsurf,'edgecolor','none');
view(2);
view(2);
title('VRAD/level 2 (m s^{-1})');
xlabel('km');
ylabel('km');

grid on
box on

applyColorScale(h1,vradPHsurf,vel_default2,colLims);

xlim(xlimits1)
ylim(ylimits1)
daspect(s1,[1 1 1]);

s2=nexttile(2);

h1=surf(XX,YY,regPHsurf,'edgecolor','none');
view(2);
view(2);
title('Regression (m s^{-1})');
xlabel('km');
ylabel('km');

grid on
box on

applyColorScale(h1,regPHsurf,vel_default2,colLims);

xlim(xlimits1)
ylim(ylimits1)
daspect(s2,[1 1 1]);

s3=nexttile(3);

h1=surf(XX,YY,regVradPHsurf,'edgecolor','none');
view(2);
view(2);
title('Regression/VRAD (m s^{-1})');
xlabel('km');
ylabel('km');

grid on
box on

applyColorScale(h1,regVradPHsurf,vel_default2,colLims);

xlim(xlimits1)
ylim(ylimits1)
daspect(s3,[1 1 1]);

linkaxes([s1,s2,s3],'xy');

print([figdir,radar,'_purpleHaze.png'],'-dpng','-r0');
=======
=======
>>>>>>> f09191c (Adding allSteps_withLegacy.m)
>>>>>>> 2323c7e (Adding allSteps_withLegacy.m)
=======
>>>>>>> 48dd00c (Adding figures for VRAD-REG paper)
