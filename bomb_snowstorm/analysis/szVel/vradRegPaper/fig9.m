% Read and diplay radar data

clear all;
close all;

addpath(genpath('~/git/lrose-test/bomb_snowstorm/analysis/'));

showPlot='on';
radar='KBGM';

figdir='/scr/cirrus1/rsfdata/projects/bomb_snowstorm/figures/vradRegPaper/';

nyquist=25.7428;

%% Infiles
regFile='/scr/cirrus1/rsfdata/projects/nexrad/tables/KBGMdoppler20181109_183814_0.48_36.08.32pt.txt';
regVradFile='/scr/cirrus1/rsfdata/projects/nexrad/matFiles/KBGM_Regression_and_VRAD_Filt_12.mat';
vradLegFile='/scr/cirrus1/rsfdata/projects/nexrad/matFiles/KBGM_VRAD_Legacy.mat';
censThreshStd=6; %7

%% Read regressino data

disp('Loading data ...')
regIn=readDataTables(regFile,' ');

%% Read mat files

regVradIn=load(regVradFile);
vradLegIn=load(vradLegFile);

%% Match azimuths

azRes=round((regIn.azimuth(2)-regIn.azimuth(1))*10)/10;
allAz=0.25:360;

ib1=[];
ib2=[];
ib3=[];
ibAll=[];
for kk=1:length(allAz)
    [minDiff1,minInd1]=min(abs(regIn.azimuth-allAz(kk)));
    [minDiff2,minInd2]=min(abs(regVradIn.az-allAz(kk)));
    [minDiff3,minInd3]=min(abs(vradLegIn.az-allAz(kk)));
    %if minDiff1<azRes/2+0.01 & minDiff2<azRes/2+0.01
        ib1=cat(1,ib1,minInd1);
        ib2=cat(1,ib2,minInd2);
        ib3=cat(1,ib3,minInd3);
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
% Level 2
% ovl=flipud(vradLegIn.ovl);
% ovl(ovl==1)=nan;
% ovl(ovl==0)=1;
vradLegIn.v1(vradLegIn.ovl==1)=nan;
lev2T=flipud((vradLegIn.v1.*vradLegIn.thr)');
lev2.VEL=nan(length(allAz),size(lev2T,2));
lev2.VEL(ibAll,:)=lev2T(ib3,:);
lev2.VEL(isinf(lev2.VEL))=nan;
% VRAD legacy
vradLegIn.vt=flipud((vradLegIn.vrad.*vradLegIn.thr)');
vradLeg.VEL=nan(length(allAz),size(vradLegIn.vt,2));
vradLeg.VEL(ibAll,:)=vradLegIn.vt(ib3,:);
vradLeg.VEL(isinf(vradLeg.VEL))=nan;

%% Cut off range
maxRange=min([regIn.range(end),regVradIn.range_km(end),vradLegIn.range_km(end)]);
minRange=max([regIn.range(1),regVradIn.range_km(1),vradLegIn.range_km(1)]);

ri1=find(reg.range>=minRange & reg.range<=maxRange);
reg.VEL=reg.VEL(:,ri1);
reg.range=reg.range(:,ri1);

ri2=find(regVradIn.range_km>=minRange & regVradIn.range_km<=maxRange);
regVrad.VEL=regVrad.VEL(:,ri2);
regInVrad.VEL=regInVrad.VEL(:,ri2);

ri3=find(vradLegIn.range_km>=minRange & vradLegIn.range_km<=maxRange);
vradLeg.VEL=vradLeg.VEL(:,ri3);
lev2.VEL=lev2.VEL(:,ri3);


%% Censor regression

regCensored=reg.VEL;

kernel=[9,5]; % Az and range of std kernel. Default: [9,5]
[stdVel,~]=fast_nd_std(reg.VEL,kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',mode(nyquist));
regCensored(stdVel>censThreshStd)=nan;

%% Fill in with regression
regVradFilled=regVrad.VEL;
regVradFilled(isnan(regVradFilled))=regCensored(isnan(regVradFilled));

%% Plot preparation

allAz(end)=359.95;
allAz(1)=0;
ang_p = deg2rad(90-allAz);

angMat=repmat(ang_p',size(reg.range,1),1);

XX1 = (reg.range.*cos(angMat));
YY1 = (reg.range.*sin(angMat));

pan1=flipud(lev2.VEL);
pan2=regVradFilled;

%% Infiles
regFile='/scr/cirrus1/rsfdata/projects/nexrad/tables/KBOX-DOPPLER20240508_131424_0.48_95.63.32PT.txt';
regVradFile='/scr/cirrus1/rsfdata/projects/nexrad/matFiles/KBOX_Regression_and_VRAD_Filt_10.mat';
vradLegFile='/scr/cirrus1/rsfdata/projects/nexrad/matFiles/KBOX_VRAD_Legacy.mat';
xlimits1=[-250,200];
ylimits1=[-200,300];

%% Read regressino data

disp('Loading data ...')
regIn=readDataTables(regFile,' ');

%% Read mat files

regVradIn=load(regVradFile);
vradLegIn=load(vradLegFile);

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
    allAz=0.25:360;
% end

ib1=[];
ib2=[];
ib3=[];
ibAll=[];
for kk=1:length(allAz)
    [minDiff1,minInd1]=min(abs(regIn.azimuth-allAz(kk)));
    [minDiff2,minInd2]=min(abs(regVradIn.az-allAz(kk)));
    [minDiff3,minInd3]=min(abs(vradLegIn.az-allAz(kk)));
    %if minDiff1<azRes/2+0.01 & minDiff2<azRes/2+0.01
        ib1=cat(1,ib1,minInd1);
        ib2=cat(1,ib2,minInd2);
        ib3=cat(1,ib3,minInd3);
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
% Level 2
% ovl=flipud(vradLegIn.ovl);
% ovl(ovl==1)=nan;
% ovl(ovl==0)=1;
vradLegIn.v1(vradLegIn.ovl==1)=nan;
lev2T=flipud((vradLegIn.v1.*vradLegIn.thr)');
lev2.VEL=nan(length(allAz),size(lev2T,2));
lev2.VEL(ibAll,:)=lev2T(ib3,:);
lev2.VEL(isinf(lev2.VEL))=nan;
% VRAD legacy
vradLegIn.vt=flipud((vradLegIn.vrad.*vradLegIn.thr)');
vradLeg.VEL=nan(length(allAz),size(vradLegIn.vt,2));
vradLeg.VEL(ibAll,:)=vradLegIn.vt(ib3,:);
vradLeg.VEL(isinf(vradLeg.VEL))=nan;

%% Cut off range
maxRange=min([regIn.range(end),regVradIn.range_km(end),vradLegIn.range_km(end)]);
minRange=max([regIn.range(1),regVradIn.range_km(1),vradLegIn.range_km(1)]);

ri1=find(reg.range>=minRange & reg.range<=maxRange);
reg.VEL=reg.VEL(:,ri1);
reg.range=reg.range(:,ri1);

ri2=find(regVradIn.range_km>=minRange & regVradIn.range_km<=maxRange);
regVrad.VEL=regVrad.VEL(:,ri2);
regInVrad.VEL=regInVrad.VEL(:,ri2);

ri3=find(vradLegIn.range_km>=minRange & vradLegIn.range_km<=maxRange);
vradLeg.VEL=vradLeg.VEL(:,ri3);
lev2.VEL=lev2.VEL(:,ri3);


%% Censor regression

censThreshStd=7; %7
regCensored=reg.VEL;

kernel=[9,5]; % Az and range of std kernel. Default: [9,5]
[stdVel,~]=fast_nd_std(reg.VEL,kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',mode(nyquist));
regCensored(stdVel>censThreshStd)=nan;

%% Fill in with regression
regVradFilled=regVrad.VEL;
regVradFilled(isnan(regVradFilled))=regCensored(isnan(regVradFilled));

%% Plot preparation
allAz(end)=359.95;
allAz(1)=0;
ang_p = deg2rad(90-allAz);

angMat=repmat(ang_p',size(reg.range,1),1);

XX2 = (reg.range.*cos(angMat));
YY2 = (reg.range.*sin(angMat));

pan3=flipud(lev2.VEL);
pan4=regVradFilled;

%% Plot

xlimits1=[-250,250];
ylimits1=[-250,250];
xlimits2=[-250,200];
ylimits2=[-200,290];

tickXY=-300:100:300;

close all
f1 = figure('Position',[200 500 930 800],'DefaultAxesFontSize',12,'visible',showPlot);
colLims=[-inf,-40,-27,-21,-17,-13,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,13,17,21,27,40,inf];

t = tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

% NEXRAD level 2
s1=nexttile(1);
hold on
h1=surf(XX1(:,1:end-1),YY1(:,1:end-1),pan1,'edgecolor','none');
view(2);
title('(a) KBGM Level-2 velocity (m s^{-1})');
%xlabel('km');
ylabel('km');

grid on
box on

applyColorScale(h1,pan1,vel_default2,colLims);

xlim(xlimits1)
ylim(ylimits1)
daspect(s1,[1 1 1]);

s1.XTick=tickXY;
s1.YTick=tickXY;

% VRAD and regression filled
s2=nexttile(2);

h1=surf(XX1,YY1,pan2,'edgecolor','none');
view(2);
title('(b) KGBM VRAD-REG filled velocity (m s^{-1})');
%xlabel('km');
%ylabel('km');

grid on
box on

applyColorScale(h1,pan2,vel_default2,colLims);

xlim(xlimits1)
ylim(ylimits1)
daspect(s2,[1 1 1]);

s2.XTick=tickXY;
s2.YTick=tickXY;

% NEXRAD level 2
s3=nexttile(3);
hold on
h1=surf(XX2(:,1:end-1),YY2(:,1:end-1),pan3,'edgecolor','none');
view(2);
title('(c) KBOX Level-2 velocity (m s^{-1})');
xlabel('km');
ylabel('km');

grid on
box on

applyColorScale(h1,pan3,vel_default2,colLims);

xlim(xlimits2)
ylim(ylimits2)
daspect(s3,[1 1 1]);

s3.XTick=tickXY;
s3.YTick=tickXY;

% VRAD and regression filled
s4=nexttile(4);

h1=surf(XX2,YY2,pan4,'edgecolor','none');
view(2);
title('(d) KBOX VRAD-REG filled velocity (m s^{-1})');
xlabel('km');
%ylabel('km');

grid on
box on

applyColorScale(h1,pan4,vel_default2,colLims);

xlim(xlimits2)
ylim(ylimits2)
daspect(s4,[1 1 1]);

s4.XTick=tickXY;
s4.YTick=tickXY;

%linkaxes([s1,s3,s4,s2,s5,s2],'xy');

print([figdir,'figure9.png'],'-dpng','-r0');

