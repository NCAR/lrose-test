% Read and diplay radar data

clear all;
close all;

% Input file
infile='/scr/cirrus1/rsfdata/projects/nexrad/tables/KDDC_LPRT_SR_20200525_030738_12pt.txt';

% Directory for figures
figdir='/scr/cirrus1/rsfdata/projects/bomb_snowstorm/figures/testFigs/';
outstr='KDDC_LPRT_SR_20200525_030738'; % First part of figure file name
% x and y limits of figures
xlimits=[-100,100];
ylimits=[-100,100];

showPlot='on'; % Display plot ('on') or just save ('off')

%% Read data
data=readDataTables(infile,' ');

%% Plot preparation

ang_p = deg2rad(90-data.azimuth);

angMat=repmat(ang_p,size(data.range,1),1);

XX = (data.range.*cos(angMat));
YY = (data.range.*sin(angMat));

%% Z
close all

figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'visible',showPlot);

s1=subplot(1,1,1);
surf(XX,YY,data.DBZ_F,'edgecolor','none');
view(2);
caxis([-10 65])
title('DBZ (dBZ)')
xlabel('km');
ylabel('km');
s1.Colormap=dbz_default2;
cb1=colorbar('XTick',-10:3:65);

grid on
box on

% Save first zoom

xlim(xlimits)
ylim(ylimits)
daspect(s1,[1 1 1]);

print([figdir,outstr,'_DBZ_zoom1.png'],'-dpng','-r0');

%% ZDR

close all

figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'visible',showPlot);

s1=subplot(1,1,1);
h=surf(XX,YY,data.ZDR_F,'edgecolor','none');
view(2);
title('ZDR (dB)')
xlabel('km');
ylabel('km');

s1.Colormap=zdr_default;
colLims=[-inf,-20,-2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.5,2,2.5,3,4,5,6,8,10,15,20,50,99,inf];
applyColorScale(h,data.ZDR_F,zdr_default,colLims);

grid on
box on

% Save first zoom

xlim(xlimits)
ylim(ylimits)
daspect(s1,[1 1 1]);

print([figdir,outstr,'_ZDR_zoom1.png'],'-dpng','-r0');

%% VEL

close all

figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'visible',showPlot);

s1=subplot(1,1,1);
h3=surf(XX,YY,data.VEL_F,'edgecolor','none');
view(2);
title('VEL (m s^{-1})')
xlabel('km');
ylabel('km');

grid on
box on

colLims=[-inf,-30,-26,-21,-17,-13,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,13,17,21,26,30,inf];
applyColorScale(h3,data.VEL_F,vel_default2,colLims);

% Save first zoom

xlim(xlimits)
ylim(ylimits)
daspect(s1,[1 1 1]);

print([figdir,outstr,'_VEL_zoom1.png'],'-dpng','-r0');

%% ORDER

if isfield(data,'REGR_ORDER')
    close all

    figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'visible',showPlot);

    s1=subplot(1,1,1);
    h=surf(XX,YY,data.REGR_ORDER,'edgecolor','none');
    view(2);
    title('ORDER')
    xlabel('km');
    ylabel('km');

    s1.Colormap=turbo(12);
    caxis([0.5,12.5]);
    colorbar;

    grid on
    box on

    % Save first zoom

    xlim(xlimits)
    ylim(ylimits)
    daspect(s1,[1 1 1]);

    print([figdir,outstr,'_ORDER_zoom1.png'],'-dpng','-r0');
end

%% PHIDP

close all

figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'visible',showPlot);

s1=subplot(1,1,1);
surf(XX,YY,data.PHIDP_F,'edgecolor','none');
view(2);
colorbar
caxis([-180,180]);
title('PHIDP (deg)')
xlabel('km');
ylabel('km');
s1.Colormap=phidp_default;

grid on
box on

% Save first zoom

xlim(xlimits)
ylim(ylimits)
daspect(s1,[1 1 1]);

print([figdir,outstr,'_PHIDP_zoom1.png'],'-dpng','-r0');

%% RHOHV

close all

figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'visible',showPlot);

s1=subplot(1,1,1);
h=surf(XX,YY,data.RHOHV_NNC_F,'edgecolor','none');
view(2);
title('RHOHV')
xlabel('km');
ylabel('km');

colLims=[-inf,0,0.7,0.8,0.85,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.975,0.98,0.985,0.99,0.995,1.1,inf];
applyColorScale(h,data.RHOHV_NNC_F,rhohv_default,colLims);

grid on
box on

% Save first zoom

xlim(xlimits)
ylim(ylimits)
daspect(s1,[1 1 1]);

print([figdir,outstr,'_RHOHV_zoom1.png'],'-dpng','-r0');

%% WIDTH

close all

figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'visible',showPlot);

s1=subplot(1,1,1);
h=surf(XX,YY,data.WIDTH_F,'edgecolor','none');
view(2);
title('WIDTH (m s^{-1})')
xlabel('km');
ylabel('km');

colLims=[-inf,0,0.5,1,1.5,2,2.5,3,4,5,6,7,8,10,12.5,15,20,25,50,inf];
applyColorScale(h,data.WIDTH_F,width_default,colLims);

grid on
box on

% Save first zoom

xlim(xlimits)
ylim(ylimits)
daspect(s1,[1 1 1]);

print([figdir,outstr,'_WIDTH_zoom1.png'],'-dpng','-r0');

%% CMD FLAG

close all

figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'visible',showPlot);

s1=subplot(1,1,1);
h=surf(XX,YY,data.CMD_FLAG,'edgecolor','none');
view(2);
title('CMD FLAG')
xlabel('km');
ylabel('km');

s1.Colormap=[0,0,1;1,0,0];
caxis([0,1]);
colorbar('Ticks',[0.25,0.75],'TickLabels',{'0','1'});

grid on
box on

% Save first zoom

xlim(xlimits)
ylim(ylimits)
daspect(s1,[1 1 1]);

print([figdir,outstr,'_CMD_zoom1.png'],'-dpng','-r0');

%% Make TRIP plot

if isfield(data,'TRIP')
    close all

    figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'visible',showPlot);

    s1=subplot(1,1,1);
    surf(XX,YY,data.TRIP,'edgecolor','none');
    view(2);
    title('TRIP')
    xlabel('km');
    ylabel('km');
    s1.Colormap=[0,0,1;1,0,0];
    caxis([0,1]);
    colorbar('Ticks',[0.25,0.75],'TickLabels',{'0','1'});

    grid on
    box on

    % Save first zoom

    xlim(xlimits)
    ylim(ylimits)
    daspect(s1,[1 1 1]);

    print([figdir,outstr,'_TRIP_zoom1.png'],'-dpng','-r0');
end