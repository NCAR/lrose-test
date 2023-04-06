% Read and diplay radar data

clear all;
close all;

addpath(genpath('~/git/lrose-test/bomb_snowstorm/analysis/utils/'));

minMaxRange=[]; % Range interval [min,max] or leave empty
minMaxAz=[]; % Azimuth interval [min,max] or leave empty
kernel=[9,5];

censorOnDBZ=1;
halfNyquist=1; % In some files the nyquist needs to be divided by 2

figdir=['/scr/cirrus1/rsfdata/projects/bomb_snowstorm/figures/statsCompare/'];

%% Loop through cases

fileID = fopen('compareFiles.txt');
inAll=textscan(fileID,'%s %s %s %f %f %f %f %f %f %f %f %s %s');
fclose(fileID);

for aa=1:size(inAll{1,1},1)

    nyquist=[];

    %% Read file 1
    infile1=inAll{1,1}(aa);

    disp(['File 1: ',infile1{:}]);

    fileType=inAll{1,12}(aa);

    if strcmp(fileType{:},'nc')
        data1in=[];

        data1in.DBZ_F=[];
        data1in.VEL_F=[];
        data1in.WIDTH_F=[];
        data1in.ZDR_F=[];
        data1in.PHIDP_F=[];
        data1in.RHOHV_NNC_F=[];
        data1in.REGR_ORDER=[];

        data1in=read_spol(infile1{:},data1in);
        nyquist=ncread(infile1{:},'nyquist_velocity');
   
    elseif strcmp(fileType{:},'table')
        data1in=readDataTables(infile1{:},' ');
        data1in.azimuth=round(data1in.azimuth);
    end

    %% Read file 2
    infile2=inAll{1,2}(aa);

    disp(['File 2: ',infile2{:}]);

    fileType=inAll{1,13}(aa);

    if strcmp(fileType{:},'nc')
        data2in=[];

        data2in.DBZ_F=[];
        data2in.VEL_F=[];
        data2in.WIDTH_F=[];
        data2in.ZDR_F=[];
        data2in.PHIDP_F=[];
        data2in.RHOHV_NNC_F=[];
        data2in.REGR_ORDER=[];

        data2in=read_spol(infile2{:},data2in);
        nyquist=ncread(infile2{:},'nyquist_velocity');

    elseif strcmp(fileType{:},'table')
        data2in=readDataTables(infile2{:},' ');
        data2in.azimuth=round(data2in.azimuth);
    end

    if isempty(nyquist)
        error('No nyquist velocity found.')
    else
        nyquist=mode(nyquist);
    end

    if halfNyquist
        nyquist=nyquist/2;
    end

    %% Match azimuths and nans

    allAz=1:360;
    data1=[];
    if ~isempty(minMaxAz)
        data1in.azimuth(data1in.azimuth<minMaxAz(1) | data1in.azimuth>minMaxAz(2))=nan;
        data2in.azimuth(data2in.azimuth<minMaxAz(1) | data2in.azimuth>minMaxAz(2))=nan;
    end
    [~,~,ib1]=intersect(allAz,data1in.azimuth);
    data1.range=data1in.range;
    data2=[];
    [~,~,ib2]=intersect(allAz,data2in.azimuth);
    data2.range=data2in.range;

    inFields=fields(data1in);
    for ii=1:size(inFields,1)
        if ~strcmp(inFields{ii},'range') & ~strcmp(inFields{ii},'time')
            data1.(inFields{ii})=nan(length(allAz),size(data1in.(inFields{ii}),2));
            data1.(inFields{ii})(ib1,:)=data1in.(inFields{ii})(ib1,:);
            data2.(inFields{ii})=nan(length(allAz),size(data2in.(inFields{ii}),2));
            data2.(inFields{ii})(ib2,:)=data2in.(inFields{ii})(ib2,:);
        end
    end

    for ii=1:size(inFields,1)
        if ~strcmp(inFields{ii},'range') & ~strcmp(inFields{ii},'time')
            % Censor on DBZ
            if censorOnDBZ & size(data1.(inFields{ii}))==size(data1.DBZ_F)
                data1.(inFields{ii})(isnan(data1.DBZ_F))=nan;
                data2.(inFields{ii})(isnan(data2.DBZ_F))=nan;
            end
            % Match nans
            data1.(inFields{ii})(isnan(data2.(inFields{ii})))=nan;
            data2.(inFields{ii})(isnan(data1.(inFields{ii})))=nan;
        end
    end

    %% Cut range
    if ~isempty(minMaxRange)

        goodInds=find(data1.range>=minMaxRange(1) & data1.range<=minMaxRange(2));

        for ii=1:size(inFields,1)
            if ~(strcmp(inFields{ii},'azimuth') | strcmp(inFields{ii},'elevation') | strcmp(inFields{ii},'time'))
                data1.(inFields{ii})=data1.(inFields{ii})(:,goodInds);
                data2.(inFields{ii})=data2.(inFields{ii})(:,goodInds);
            end
        end
    end

    %% Plot preparation

    ang_p = deg2rad(90-data1.azimuth);

    angMat=repmat(ang_p,size(data1.range,1),1);

    XX = (data1.range.*cos(angMat));
    YY = (data1.range.*sin(angMat));

    %% Loop through fields
    for jj=5:10
        %% Standard deviations

        if strcmp(inFields{jj},'VEL_F')
            [stdVar1,~]=fast_nd_std(data1.(inFields{jj}),kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',nyquist);
            [stdVar2,~]=fast_nd_std(data2.(inFields{jj}),kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',nyquist);
        else
            [stdVar1,~]=fast_nd_std(data1.(inFields{jj}),kernel,'mode','partial','nan_std',1);
            [stdVar2,~]=fast_nd_std(data2.(inFields{jj}),kernel,'mode','partial','nan_std',1);
        end

        stdVar1(isnan(data1.(inFields{jj})))=nan;
        stdVar2(isnan(data2.(inFields{jj})))=nan;

        %% Plot
        close all

        figure('Position',[200 500 2200 1200],'DefaultAxesFontSize',12);
        colormap('jet');

        s1=subplot(2,3,1);
        pBottom=prctile(data1.(inFields{jj}),1,'all');
        pTop=prctile(data1.(inFields{jj}),99,'all');
        surf(XX,YY,data1.(inFields{jj}),'edgecolor','none');
        view(2);
        caxis([pBottom,pTop]);
        colorbar;
        title([inFields{jj},' file 1'],'Interpreter','none');
        xlabel('km');
        ylabel('km');

        grid on
        box on

        s2=subplot(2,3,2);
        surf(XX,YY,data2.(inFields{jj}),'edgecolor','none');
        view(2);
        caxis([pBottom,pTop]);
        colorbar;
        title([inFields{jj},' file 2'],'Interpreter','none');
        xlabel('km');
        ylabel('km');

        grid on
        box on

        s3=subplot(2,3,3);
        diffField=data2.(inFields{jj})-data1.(inFields{jj});
        pBottom=prctile(diffField,1,'all');
        pTop=prctile(diffField,99,'all');
        lim=max(abs([pTop,pBottom]));
        surf(XX,YY,diffField,'edgecolor','none');
        view(2);
        caxis([-lim,lim]);
        s3.Colormap=velCols;
        colorbar;
        title([inFields{jj},' file 2 - file 1'],'Interpreter','none');
        xlabel('km');
        ylabel('km');
        
        grid on
        box on

        s4=subplot(2,3,4);
        pBottom=prctile(stdVar1,10,'all');
        pTop=prctile(stdVar1,90,'all');
        surf(XX,YY,stdVar1,'edgecolor','none');
        view(2);
        caxis([pBottom,pTop]);
        colorbar;
        title(['std file 1'],'Interpreter','none');
        xlabel('km');
        ylabel('km');
        
        grid on
        box on

        s5=subplot(2,3,5);
        surf(XX,YY,stdVar2,'edgecolor','none');
        view(2);
        caxis([pBottom,pTop]);
        colorbar;
        title(['std file 2'],'Interpreter','none');
        xlabel('km');
        ylabel('km');
        
        grid on
        box on

        s6=subplot(2,3,6);
        diffField=stdVar2-stdVar1;
        pBottom=prctile(diffField,1,'all');
        pTop=prctile(diffField,99,'all');
        lim=max(abs([pTop,pBottom]));
        surf(XX,YY,diffField,'edgecolor','none');
        view(2);
        caxis([-lim,lim]);
        s6.Colormap=velCols;
        colorbar;
        title(['std file 2 - std file 1'],'Interpreter','none');
        xlabel('km');
        ylabel('km');
        
        grid on
        box on

        outstr=inAll{1,3}(aa);
        outstr=outstr{:};
        mtit([outstr],'fontsize',14,'xoff',0,'yoff',0.05,'interpreter','none');

        %% Save first zoom

        if ~isnan(minMaxRange)
            outstr=[outstr,'_range',num2str(minMaxRange(1)),'to',num2str(minMaxRange(2))];
        end
        if ~isnan(minMaxAz)
            outstr=[outstr,'_az',num2str(minMaxAz(1)),'to',num2str(minMaxAz(2))];
        end

        mkdir(figdir,outstr);

        linkaxes([s1,s2,s3,s4,s5,s6],'xy');        

        xlimits1=[inAll{1,4}(aa),inAll{1,5}(aa)];
        ylimits1=[inAll{1,6}(aa),inAll{1,7}(aa)];

        xlim(xlimits1)
        ylim(ylimits1)
        daspect(s1,[1 1 1]);
        daspect(s2,[1 1 1]);
        daspect(s3,[1 1 1]);
        daspect(s4,[1 1 1]);
        daspect(s5,[1 1 1]);
        daspect(s6,[1 1 1]);
       
        print([figdir,outstr,'/',outstr,'_',inFields{jj},'_zoom1.png'],'-dpng','-r0');

        %% Save second zoom

        xlimits2=[inAll{1,8}(aa),inAll{1,9}(aa)];
        ylimits2=[inAll{1,10}(aa),inAll{1,11}(aa)];

        xlim(xlimits2)
        ylim(ylimits2)
        daspect(s1,[1 1 1]);
        daspect(s2,[1 1 1]);
        daspect(s3,[1 1 1]);
        daspect(s4,[1 1 1]);
        daspect(s5,[1 1 1]);
        daspect(s6,[1 1 1]);
       
        print([figdir,outstr,'/',outstr,'_',inFields{jj},'_zoom2.png'],'-dpng','-r0');
    end
end