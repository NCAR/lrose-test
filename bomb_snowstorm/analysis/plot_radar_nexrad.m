% Read and diplay radar data

clear all;
close all;

addpath(genpath('~/git/lrose-test/bomb_snowstorm/analysis/utils/'));

maxRange=[];

showPlot='on';
saveData=1;

%% Loop through cases

fileID = fopen('plotFiles_nexrad.txt');
inAll=textscan(fileID,'%s %s %f %f %f %f %f %f %f %f %s %s %f');
fclose(fileID);

for aa=2:size(inAll{1,1},1)

    infile=inAll{1,1}(aa);

    disp(['File ',num2str(aa), ' of ',num2str(size(inAll{1,1},1))]);
    disp(infile{:});

    inst=inAll{1,12}(aa);
    if strcmp(inst{:},'bs')
        figdir=['/scr/cirrus1/rsfdata/projects/bomb_snowstorm/figures/'];
    elseif strcmp(inst{:},'kddc')
        figdir=['/scr/cirrus1/rsfdata/projects/nexrad/figures/kddc/'];
    elseif strcmp(inst{:},'kftg')
        figdir=['/scr/cirrus1/rsfdata/projects/nexrad/figures/kftg/'];
    end

    fileType=inAll{1,11}(aa);

    data=[];

    data.REF=[];
    data.VEL=[];
    data.SW=[];
    data.ZDR=[];
    data.PHI=[];
    data.RHO=[];
    
    data=read_spol(infile{:},data);

    data=data(inAll{1,13}(aa));
        
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

    %% Plot preparation

    ang_p = deg2rad(90-data.azimuth);

    angMat=repmat(ang_p,size(data.range,1),1);

    XX = (data.range.*cos(angMat));
    YY = (data.range.*sin(angMat));

    %% Z
    close all

    figure('Position',[200 500 2100 1200],'DefaultAxesFontSize',12,'visible',showPlot);

    s1=subplot(2,3,1);
    surf(XX,YY,data.REF,'edgecolor','none');
    view(2);
    caxis([-10 65])
    title('DBZ (dBZ)')
    xlabel('km');
    ylabel('km');
    s1.Colormap=dbz_default2;
    cb1=colorbar('XTick',-10:3:65);

    grid on
    box on
                
    %% ZDR

    s2=subplot(2,3,2);
    h=surf(XX,YY,data.ZDR,'edgecolor','none');
    view(2);
    title('ZDR (dB)')
    xlabel('km');
    ylabel('km');

    s2.Colormap=zdr_default;
    colLims=[-inf,-20,-2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.5,2,2.5,3,4,5,6,8,10,15,20,50,99,inf];
    applyColorScale(h,data.ZDR,zdr_default,colLims);

    grid on
    box on

    %% VEL

    s3=subplot(2,3,3);
    h3=surf(XX,YY,data.VEL,'edgecolor','none');
    view(2);
    title('VEL (m s^{-1})')
    xlabel('km');
    ylabel('km');

    grid on
    box on

    colLims=[-inf,-30,-26,-21,-17,-13,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,13,17,21,26,30,inf];
    applyColorScale(h3,data.VEL,vel_default2,colLims);

    %% PHIDP

    if strcmp(inst{:},'kddc') & strcmp(fileType{:},'nc')
        data.PHI=wrapTo360(data.PHI);
        data.PHI=data.PHI-90;
    end

    s4=subplot(2,3,4);
    surf(XX,YY,data.PHI,'edgecolor','none');
    view(2);
    colorbar
    if strcmp(inst{:},'bs')
        caxis([-60,92]);
    else
        caxis([-180,180]);
    end
    title('PHIDP (deg)')
    xlabel('km');
    ylabel('km');
    s4.Colormap=phidp_default;

    grid on
    box on

    %% RHOHV

    s5=subplot(2,3,5);
    h=surf(XX,YY,data.RHO,'edgecolor','none');
    view(2);
    title('RHOHV')
    xlabel('km');
    ylabel('km');

    colLims=[-inf,0,0.7,0.8,0.85,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.975,0.98,0.985,0.99,0.995,1.1,inf];
    applyColorScale(h,data.RHO,rhohv_default,colLims);

    grid on
    box on

    %% WIDTH

    s6=subplot(2,3,6);
    h=surf(XX,YY,data.SW,'edgecolor','none');
    view(2);
    title('WIDTH (m s^{-1})')
    xlabel('km');
    ylabel('km');

    colLims=[-inf,0,0.5,1,1.5,2,2.5,3,4,5,6,7,8,10,12.5,15,20,25,50,inf];
    applyColorScale(h,data.SW,width_default,colLims);

    grid on
    box on

    %% Save first zoom

    linkaxes([s1,s2,s3,s4,s5,s6],'xy');
    outstr=inAll{1,2}(aa);
    outstr=outstr{:};

    xlimits1=[inAll{1,3}(aa),inAll{1,4}(aa)];
    ylimits1=[inAll{1,5}(aa),inAll{1,6}(aa)];
    
    xlim(xlimits1)
    ylim(ylimits1)
    daspect(s1,[1 1 1]);
    daspect(s2,[1 1 1]);
    daspect(s3,[1 1 1]);
    daspect(s4,[1 1 1]);
    daspect(s5,[1 1 1]);
    daspect(s6,[1 1 1]);
        
    print([figdir,outstr,'_zoom1.png'],'-dpng','-r0');

    %% Save second zoom

    xlimits2=[inAll{1,7}(aa),inAll{1,8}(aa)];
    ylimits2=[inAll{1,9}(aa),inAll{1,10}(aa)];
    
    xlim(xlimits2)
    ylim(ylimits2)
    daspect(s1,[1 1 1]);
    daspect(s2,[1 1 1]);
    daspect(s3,[1 1 1]);
    daspect(s4,[1 1 1]);
    daspect(s5,[1 1 1]);
    daspect(s6,[1 1 1]);
           
    print([figdir,outstr,'_zoom2.png'],'-dpng','-r0');

    if saveData
        save([figdir,outstr,'.mat'],'data');
    end

end