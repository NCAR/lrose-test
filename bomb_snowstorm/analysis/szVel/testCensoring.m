% Read and diplay radar data

clear all;
close all;

addpath(genpath('~/git/lrose-test/bomb_snowstorm/analysis/'));

showPlot='on';
figdir='/scr/cirrus1/rsfdata/projects/nexrad/figures/szComp/censoring/';

%% Infiles
fileList=readtable('regFileList.txt','Delimiter',' ');

nyquist=28.39;

%% Read Data
for ii=1:size(fileList,1)

    disp(['Loading file ',num2str(ii),' of ',num2str(size(fileList,1)),'...']);

    regFile=fileList.File{ii};
    dataReg=readDataTables(regFile,' ');

    %% Censor regression

    kernel=[9,5]; % Az and range of std kernel. Default: [9,5]
    [stdVel,~]=fast_nd_std(dataReg.VEL_F,kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',mode(nyquist));

    regVelCensored=dataReg.VEL_F;
    regVelCensored(stdVel>9)=nan;

    %% Plot preparation

    ang_p = deg2rad(90-dataReg.azimuth);

    angMat=repmat(ang_p,size(dataReg.range,1),1);

    XX = (dataReg.range.*cos(angMat));
    YY = (dataReg.range.*sin(angMat));

    xlimits1=[-250,250];
    ylimits1=[-250,250];

    %% Plot

    close all
    f1 = figure('Position',[200 500 1200 500],'DefaultAxesFontSize',12,'visible',showPlot);

    t = tiledlayout(1,2,'TileSpacing','tight','Padding','tight');

    s1=nexttile(1);

    h1=surf(XX,YY,dataReg.VEL_F,'edgecolor','none');
    view(2);
    title('VEL regression (m s^{-1})');
    xlabel('km');
    ylabel('km');

    grid on
    box on

    colLims=[-inf,-30,-26,-21,-17,-13,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,13,17,21,26,30,inf];
    applyColorScale(h1,dataReg.VEL_F,vel_default2,colLims);

    xlim(xlimits1)
    ylim(ylimits1)
    daspect(s1,[1 1 1]);

    s2=nexttile(2);

    h1=surf(XX,YY,regVelCensored,'edgecolor','none');
    view(2);
    title('VEL regression censored (m s^{-1})');
    xlabel('km');
    ylabel('km');

    grid on
    box on

    colLims=[-inf,-30,-26,-21,-17,-13,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,13,17,21,26,30,inf];
    applyColorScale(h1,regVelCensored,vel_default2,colLims);

    xlim(xlimits1)
    ylim(ylimits1)
    daspect(s2,[1 1 1]);

    regSpl=strsplit(regFile,'/');
    outName=regSpl{8}(1:end-3);
    print([figdir,outName,'png'],'-dpng','-r0');
end