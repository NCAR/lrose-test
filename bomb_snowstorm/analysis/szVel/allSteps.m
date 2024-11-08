% Read and diplay radar data

clear all;
close all;

addpath(genpath('~/git/lrose-test/bomb_snowstorm/analysis/'));

showPlot='on';
radar='KFTG';

figdir=['/scr/cirrus1/rsfdata/projects/nexrad/figures/szComp/',radar,'/'];

nyquist=28.39;

%% Infiles
if strcmp(radar,'KFTG')
    regFile='/scr/cirrus1/rsfdata/projects/nexrad/tables/KFTG_SZ_20220329_190532_0.48_272.92_O37_64pts_V5.txt';
    regVradFile='/scr/cirrus1/rsfdata/projects/nexrad/matFiles/KFTG_Regression_and_VRAD.mat';
    vradLegFile='/scr/cirrus1/rsfdata/projects/nexrad/matFiles/KFTG_VRAD_Legacy.mat';
elseif strcmp(radar,'KLWX')
end

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
    allAz=1:360;
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

reg.VEL=nan(length(allAz),size(regIn.VEL_F,2));
reg.VEL(ibAll,:)=regIn.VEL_F(ib1,:);
regVradIn.vt=(regVradIn.v1.*regVradIn.thr)';
regVrad.VEL=nan(length(allAz),size(regVradIn.vt,2));
regVrad.VEL(ibAll,:)=regVradIn.vt(ib2,:);
vradLegIn.vt=(vradLegIn.v1.*vradLegIn.thr)';
vradLeg.VEL=nan(length(allAz),size(vradLegIn.vt,2));
vradLeg.VEL(ibAll,:)=vradLegIn.vt(ib3,:);

%% Cut off range
maxRange=min([regIn.range(end),regVradIn.range_km(end),vradLegIn.range_km(end)]);
minRange=max([regIn.range(1),regVradIn.range_km(1),vradLegIn.range_km(1)]);

ri1=find(reg.range>=minRange & reg.range<=maxRange);
reg.VEL=reg.VEL(:,ri1);
reg.range=reg.range(:,ri1);

ri2=find(regVradIn.range_km>=minRange & regVradIn.range_km<=maxRange);
regVrad.VEL=regVrad.VEL(:,ri2);

ri3=find(vradLegIn.range_km>=minRange & vradLegIn.range_km<=maxRange);
vradLeg.VEL=vradLeg.VEL(:,ri3);


%% Censor regression

kernel=[9,5]; % Az and range of std kernel. Default: [9,5]
[stdVel,~]=fast_nd_std(reg.VEL,kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',mode(nyquist));

regCensored=reg.VEL;
regCensored(stdVel>9)=nan;

%% Fill in with regression
regVradFilled=regVrad.VEL;
regVradFilled(isnan(regVradFilled))=regCensored(isnan(regVradFilled));

%% Plot preparation

xlimits1=[-200,260];
ylimits1=[-220,220];

ang_p = deg2rad(90-allAz);

angMat=repmat(ang_p',size(reg.range,1),1);

XX = (reg.range.*cos(angMat));
YY = (reg.range.*sin(angMat));

%% Plot

close all
f1 = figure('Position',[200 500 1800 1000],'DefaultAxesFontSize',12,'visible',showPlot);

t = tiledlayout(2,3,'TileSpacing','tight','Padding','tight');

% Regression original
s2=nexttile(2);

h1=surf(XX,YY,reg.VEL,'edgecolor','none');
view(2);
title('VEL regression (m s^{-1})');
xlabel('km');
ylabel('km');

grid on
box on

colLims=[-inf,-30,-26,-21,-17,-13,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,13,17,21,26,30,inf];
applyColorScale(h1,reg.VEL,vel_default2,colLims);

xlim(xlimits1)
ylim(ylimits1)
daspect(s2,[1 1 1]);

% Regression censored
s3=nexttile(3);

h1=surf(XX,YY,regCensored,'edgecolor','none');
view(2);
title('VEL regression censored (m s^{-1})');
xlabel('km');
ylabel('km');

grid on
box on

colLims=[-inf,-30,-26,-21,-17,-13,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,13,17,21,26,30,inf];
applyColorScale(h1,regCensored,vel_default2,colLims);

xlim(xlimits1)
ylim(ylimits1)
daspect(s3,[1 1 1]);

% VRAD legacy
s4=nexttile(4);

h1=surf(XX,YY,flipud(vradLeg.VEL),'edgecolor','none');
view(2);
title('VEL VRAD legacy (m s^{-1})');
xlabel('km');
ylabel('km');

grid on
box on

colLims=[-inf,-30,-26,-21,-17,-13,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,13,17,21,26,30,inf];
applyColorScale(h1,flipud(vradLeg.VEL),vel_default2,colLims);

xlim(xlimits1)
ylim(ylimits1)
daspect(s4,[1 1 1]);

% VRAD and regression
s5=nexttile(5);

h1=surf(XX,YY,regVrad.VEL,'edgecolor','none');
view(2);
title('VEL VRAD and regression legacy (m s^{-1})');
xlabel('km');
ylabel('km');

grid on
box on

colLims=[-inf,-30,-26,-21,-17,-13,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,13,17,21,26,30,inf];
applyColorScale(h1,regVrad.VEL,vel_default2,colLims);

xlim(xlimits1)
ylim(ylimits1)
daspect(s5,[1 1 1]);

% VRAD and regression filled
s6=nexttile(6);

h1=surf(XX,YY,regVradFilled,'edgecolor','none');
view(2);
title('VEL VRAD and regression legacy (m s^{-1})');
xlabel('km');
ylabel('km');

grid on
box on

colLims=[-inf,-30,-26,-21,-17,-13,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,13,17,21,26,30,inf];
applyColorScale(h1,regVradFilled,vel_default2,colLims);

xlim(xlimits1)
ylim(ylimits1)
daspect(s6,[1 1 1]);

print([figdir,'VEL_allSteps.png'],'-dpng','-r0');