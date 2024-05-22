% Read and diplay radar data

clear all;
close all;

addpath(genpath('~/git/lrose-test/bomb_snowstorm/analysis/utils/'));

maxRange=[];

showPlot='on';

%% Loop through cases

fileID = fopen('plotFiles_nexrad_indVars.txt');
inAll=textscan(fileID,'%s %s %f %f %f %f %f %f %f %f %s %s %f');
fclose(fileID);

for aa=8:size(inAll{1,1},1)

    infile=inAll{1,1}(aa);

    disp(['File ',num2str(aa), ' of ',num2str(size(inAll{1,1},1))]);
    disp(infile{:});

    inst=inAll{1,12}(aa);
    if strcmp(inst{:},'bs')
        figdirBase=['/scr/cirrus1/rsfdata/projects/bomb_snowstorm/figures/'];
    elseif strcmp(inst{:},'kddc')
        figdirBase=['/scr/cirrus1/rsfdata/projects/nexrad/figures/kddc/'];
    elseif strcmp(inst{:},'kftg')
        figdirBase=['/scr/cirrus1/rsfdata/projects/nexrad/figures/kftg/'];
    end

    fileType=inAll{1,11}(aa);

    data=[];

    data.REF=[];
    data.DBZ=[];
    data.VEL=[];
    data.SW=[];
    data.WIDTH=[];
    data.ZDR=[];
    data.PHI=[];
    data.PHIDP=[];
    data.RHO=[];
    data.RHOHV=[];
    data.PURPLE_HAZE=[];
    
    data=read_spol(infile{:},data);

    data=data(inAll{1,13}(aa));

    if isfield(data,'REF')
        data.DBZ=data.REF;
    end
    if isfield(data,'SW')
        data.WIDTH=data.SW;
    end
    if isfield(data,'PHI')
        data.PHIDP=data.PHI;
    end
    if isfield(data,'RHO')
        data.RHOHV=data.RHO;
    end

    %% Cut range
    if ~isempty(maxRange)
        inFields=fields(data);

        goodInds=find(data.range<=maxRange);

        for ii=1:size(inFields,1)
            if ~(strcmp(inFields{ii},'azimuth') | strcmp(inFields{ii},'elevation') | strcmp(inFields{ii},'time'))
                data.(inFields{ii})=data.(inFields{ii})(:,goodInds);
            end
        end
    end

%% Add purple haze if exists

    % if isfield(data,'PURPLE_HAZE')
    %     inFields=fields(data);
    %     for ii=1:size(inFields,1)
    %         if ~(strcmp(inFields{ii},'azimuth') | strcmp(inFields{ii},'elevation') | strcmp(inFields{ii},'time') ...
    %                 | strcmp(inFields{ii},'range') | strcmp(inFields{ii},'PURPLE_HAZE'))
    %             data.(inFields{ii})(data.PURPLE_HAZE==1)=99999;
    %         end
    %     end
    % end
    %% Plot preparation

    ang_p = deg2rad(90-data.azimuth);

    angMat=repmat(ang_p,size(data.range,1),1);

    XX = (data.range.*cos(angMat));
    YY = (data.range.*sin(angMat));

    outstr=inAll{1,2}(aa);
    outstr=outstr{:};

    figdir=[figdirBase,outstr,'/'];
    mkdir(figdir);

    xlimits1=[inAll{1,3}(aa),inAll{1,4}(aa)];
    ylimits1=[inAll{1,5}(aa),inAll{1,6}(aa)];

    xlimits2=[inAll{1,7}(aa),inAll{1,8}(aa)];
    ylimits2=[inAll{1,9}(aa),inAll{1,10}(aa)];

  
    %% Z
    close all

    figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'visible',showPlot);

    s1=subplot(1,1,1);
    surf(XX,YY,data.DBZ,'edgecolor','none');
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

    xlim(xlimits1)
    ylim(ylimits1)
    daspect(s1,[1 1 1]);

    print([figdir,outstr,'_DBZ_zoom1.png'],'-dpng','-r0');

    % Save second zoom

    xlim(xlimits2)
    ylim(ylimits2)
    daspect(s1,[1 1 1]);

    print([figdir,outstr,'_DBZ_zoom2.png'],'-dpng','-r0');

    %% ZDR

    close all

    figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'visible',showPlot);

    s1=subplot(1,1,1);
    h=surf(XX,YY,data.ZDR,'edgecolor','none');
    view(2);
    title('ZDR (dB)')
    xlabel('km');
    ylabel('km');

    s1.Colormap=zdr_default;
    colLims=[-inf,-20,-2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.5,2,2.5,3,4,5,6,8,10,15,20,50,99,inf];
    applyColorScale(h,data.ZDR,zdr_default,colLims);

    grid on
    box on

    % Save first zoom

    xlim(xlimits1)
    ylim(ylimits1)
    daspect(s1,[1 1 1]);

    print([figdir,outstr,'_ZDR_zoom1.png'],'-dpng','-r0');

    % Save second zoom

    xlim(xlimits2)
    ylim(ylimits2)
    daspect(s1,[1 1 1]);

    print([figdir,outstr,'_ZDR_zoom2.png'],'-dpng','-r0');

    %% VEL

    close all

    figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'visible',showPlot);

    s1=subplot(1,1,1);
    h3=surf(XX,YY,data.VEL,'edgecolor','none');
    view(2);
    title('VEL (m s^{-1})')
    xlabel('km');
    ylabel('km');

    grid on
    box on

    colLims=[-inf,-30,-26,-21,-17,-13,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,13,17,21,26,30,inf];
    applyColorScale(h3,data.VEL,vel_default2,colLims);

    % Save first zoom

    xlim(xlimits1)
    ylim(ylimits1)
    daspect(s1,[1 1 1]);

    print([figdir,outstr,'_VEL_zoom1.png'],'-dpng','-r0');

    % Save second zoom

    xlim(xlimits2)
    ylim(ylimits2)
    daspect(s1,[1 1 1]);

    print([figdir,outstr,'_VEL_zoom2.png'],'-dpng','-r0');

    %% PHIDP

    if strcmp(inst{:},'kddc') & strcmp(fileType{:},'nc')
        data.PHIDP=wrapTo360(data.PHIDP);
        data.PHIDP=data.PHIDP-90;
    end

    close all

    figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'visible',showPlot);

    s1=subplot(1,1,1);
    surf(XX,YY,data.PHIDP,'edgecolor','none');
    view(2);
    colorbar
    if strcmp(inst{:},'bs')
        caxis([-60,92]);
    elseif strcmp(inst{:},'kftg')
        caxis([0,114]);
    else
        caxis([-180,180]);
    end
    title('PHIDP (deg)')
    xlabel('km');
    ylabel('km');
    s1.Colormap=phidp_default;

    grid on
    box on

    % Save first zoom

    xlim(xlimits1)
    ylim(ylimits1)
    daspect(s1,[1 1 1]);

    print([figdir,outstr,'_PHIDP_zoom1.png'],'-dpng','-r0');

    % Save second zoom

    xlim(xlimits2)
    ylim(ylimits2)
    daspect(s1,[1 1 1]);

    print([figdir,outstr,'_PHIDP_zoom2.png'],'-dpng','-r0');

    %% RHOHV

    close all

    figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'visible',showPlot);

    s1=subplot(1,1,1);
    h=surf(XX,YY,data.RHOHV,'edgecolor','none');
    view(2);
    title('RHOHV')
    xlabel('km');
    ylabel('km');

    colLims=[-inf,0,0.7,0.8,0.85,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.975,0.98,0.985,0.99,0.995,1.1,inf];
    applyColorScale(h,data.RHOHV,rhohv_default,colLims);

    grid on
    box on

    % Save first zoom

    xlim(xlimits1)
    ylim(ylimits1)
    daspect(s1,[1 1 1]);

    print([figdir,outstr,'_RHOHV_zoom1.png'],'-dpng','-r0');

    % Save second zoom

    xlim(xlimits2)
    ylim(ylimits2)
    daspect(s1,[1 1 1]);

    print([figdir,outstr,'_RHOHV_zoom2.png'],'-dpng','-r0');

    %% WIDTH

    close all

    figure('Position',[200 500 1000 800],'DefaultAxesFontSize',12,'visible',showPlot);

    s1=subplot(1,1,1);
    h=surf(XX,YY,data.WIDTH,'edgecolor','none');
    view(2);
    title('WIDTH (m s^{-1})')
    xlabel('km');
    ylabel('km');

    colLims=[-inf,0,0.5,1,1.5,2,2.5,3,4,5,6,7,8,10,12.5,15,20,25,50,inf];
    applyColorScale(h,data.WIDTH,width_default,colLims);

    grid on
    box on

    % Save first zoom

    xlim(xlimits1)
    ylim(ylimits1)
    daspect(s1,[1 1 1]);

    print([figdir,outstr,'_WIDTH_zoom1.png'],'-dpng','-r0');

    % Save second zoom

    xlim(xlimits2)
    ylim(ylimits2)
    daspect(s1,[1 1 1]);

    print([figdir,outstr,'_WIDTH_zoom2.png'],'-dpng','-r0');

end