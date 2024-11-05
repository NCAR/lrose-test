% Read and diplay radar data

clear all;
close all;

addpath(genpath('~/git/lrose-test/bomb_snowstorm/analysis/utils/'));

minMaxRangeOrig=[]; % Range interval [min,max] or leave empty
minMaxAz=[]; % Azimuth interval [min,max] or leave empty
kernel=[9,5]; % Az and range of std kernel. Default: [9,5]

censorOnDBZ=0;
censorOnVEL=0;
censorOnCMD=1;
censorOnTRIP=0; % Only use weak trip (0).
%%%%%%%%%%%%%%
tripToSnr=0; % The last (10th) variable that is read in John's files is TRIP. Sometimes it is actually SNR.
censorOnSNR=[]; % Set to empty if not used !!!!!!! Only use areas with SNR above XX dB
%%%%%%%%%%%%%%
halfNyquist=0; % In some files the nyquist needs to be divided by 2
removeZeros=0;

%% Loop through cases

fileID = fopen('compareFiles_cases2425.txt');
inAll=textscan(fileID,'%s %s %s %f %f %f %f %f %f %f %f %s %s %s %f');
fclose(fileID);

showPlot='on';

for aa=1:size(inAll{1,1},1)

    nyquist=[];

    %% Read file 1
    infile1=inAll{1,1}(aa);

    disp(['File 1: ',infile1{:}]);

    figdir=['/scr/cirrus1/rsfdata/projects/nexrad/figures/cases2425/hist2D/'];    

    fileType=inAll{1,12}(aa);

    if strcmp(fileType{:},'nc')
        data1in=[];

        data1in.DBZ_F=[];
        data1in.VEL_F=[];
        data1in.WIDTH_F=[];
        data1in.ZDR_F=[];
        data1in.PHIDP_F=[];
        data1in.RHOHV_F=[];
        data1in.REGR_FILT_POLY_ORDER=[];
        data1in.CMD_FLAG=[];

        data1in=read_spol(infile1{:},data1in);
        nyquist=ncread(infile1{:},'nyquist_velocity');

    elseif strcmp(fileType{:},'nexrad')
        data1in=[];

        data1in.DBZ=[];
        data1in.VEL=[];
        data1in.WIDTH=[];
        data1in.ZDR=[];
        data1in.PHIDP=[];
        data1in.RHOHV=[];

        data1in=read_spol(infile1{:},data1in);
        nyquist=ncread(infile1{:},'nyquist_velocity');

        data1in=data1in(inAll{1,15}(aa));

        data1in.DBZ_F=data1in.DBZ;
        data1in.VEL_F=data1in.VEL;
        data1in.WIDTH_F=data1in.WIDTH;
        data1in.ZDR_F=data1in.ZDR;
        data1in.PHIDP_F=data1in.PHIDP;
        data1in.RHOHV_F=data1in.RHOHV;

    elseif strcmp(fileType{:},'nexradLevel2')
        data1in=[];

        data1in.REF=[];
        data1in.VEL=[];
        data1in.SW=[];
        data1in.ZDR=[];
        data1in.PHI=[];
        data1in.RHO=[];

        data1in=read_spol(infile1{:},data1in);
        nyquist=ncread(infile1{:},'nyquist_velocity');

        data1in=data1in(inAll{1,15}(aa));

        data1in.DBZ_F=data1in.REF;
        data1in.VEL_F=data1in.VEL;
        data1in.WIDTH_F=data1in.SW;
        data1in.ZDR_F=data1in.ZDR;
        data1in.PHIDP_F=data1in.PHI;
        data1in.RHOHV_F=data1in.RHO;

    elseif strcmp(fileType{:},'table')
        data1in=readDataTables(infile1{:},' ');
        data1in.azimuth=round(data1in.azimuth);
        data1in.RHOHV_F=data1in.RHOHV_NNC_F;
        if tripToSnr
            data1in.SNR=data1in.TRIP;
            data1in=rmfield(data1in,'TRIP');
        end
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
        data2in.RHOHV_F=[];
        data2in.REGR_FILT_POLY_ORDER=[];
        data2in.CMD_FLAG=[];

        data2in=read_spol(infile2{:},data2in);
        nyquist=ncread(infile2{:},'nyquist_velocity');

    elseif strcmp(fileType{:},'nexrad')
        data2in=[];

        data2in.DBZ=[];
        data2in.VEL=[];
        data2in.WIDTH=[];
        data2in.ZDR=[];
        data2in.PHIDP=[];
        data2in.RHOHV=[];

        data1in=read_spol(infile2{:},data2in);
        nyquist=ncread(infile2{:},'nyquist_velocity');

        data2in=data2in(inAll{1,15}(aa));

        data2in.DBZ_F=data2in.DBZ;
        data2in.VEL_F=data2in.VEL;
        data2in.WIDTH_F=data2in.WIDTH;
        data2in.ZDR_F=data2in.ZDR;
        data2in.PHIDP_F=data2in.PHIDP;
        data2in.RHOHV_F=data2in.RHOHV;

    elseif strcmp(fileType{:},'nexradLevel2')
        data2in=[];

        data2in.REF=[];
        data2in.VEL=[];
        data2in.SW=[];
        data2in.ZDR=[];
        data2in.PHI=[];
        data2in.RHO=[];

        data2in=read_spol(infile1{:},data2in);
        nyquist=ncread(infile1{:},'nyquist_velocity');

        data2in=data2in(inAll{1,15}(aa));

        data2in.DBZ_F=data2in.REF;
        data2in.VEL_F=data2in.VEL;
        data2in.WIDTH_F=data2in.SW;
        data2in.ZDR_F=data2in.ZDR;
        data2in.PHIDP_F=data2in.PHI;
        data2in.RHOHV_F=data2in.RHO;

    elseif strcmp(fileType{:},'table')
        data2in=readDataTables(infile2{:},' ');
        %data2in.azimuth=round(data2in.azimuth);
        if tripToSnr
            data2in.SNR=data2in.TRIP;
            data2in=rmfield(data2in,'TRIP');
        end
        data2in.RHOHV_F=data2in.RHOHV_NNC_F;
    elseif strcmp(fileType{:},'mat')
        load(infile2{:});
        addnan=nan(size(data.REF,1),8);
        data2in.DBZ_F=cat(2,addnan,data.REF(:,1:end-1));
        data2in.VEL_F=cat(2,addnan,data.VEL(:,1:end-1));
        data2in.WIDTH_F=cat(2,addnan,data.SW(:,1:end-1));
        data2in.ZDR_F=cat(2,addnan,data.ZDR(:,1:end-1));
        data2in.PHIDP_F=cat(2,addnan,data.PHI(:,1:end-1));
        data2in.azimuth=data.azimuth;
        data2in.range=data1in.range;

        data2in.azimuth=floor(data2in.azimuth);
        data2in.azimuth(2:2:length(data2in.azimuth))=data2in.azimuth(2:2:length(data2in.azimuth))+0.5;
    end

    if isempty(nyquist)
        warning('No nyquist velocity found. Using 26.675 m/s')
        nyquist=26.675;
    else
        nyquist=mode(nyquist);
    end

    if halfNyquist
        nyquist=nyquist/2;
    end

    %% Cut range
    inFields1=fields(data1in);
    inFields2=fields(data2in);
    inFields=intersect(inFields1,inFields2);

    if isempty(minMaxRangeOrig)
        minMaxRange=max([data1in.range(1),data2in.range(1)]);
        minMaxRange=[minMaxRange,min([data1in.range(end),data2in.range(end)])];
    else
        minMaxRange=minMaxRangeOrig;
    end
    goodInds1=find(data1in.range>=minMaxRange(1)-0.001 & data1in.range<=minMaxRange(2)+0.001);
    goodInds2=find(data2in.range>=minMaxRange(1)-0.001 & data2in.range<=minMaxRange(2)+0.001);

    for ii=1:size(inFields,1)
        if ~(strcmp(inFields{ii},'azimuth') | strcmp(inFields{ii},'elevation') | strcmp(inFields{ii},'time'))
            data1in.(inFields{ii})=data1in.(inFields{ii})(:,goodInds1);
            data2in.(inFields{ii})=data2in.(inFields{ii})(:,goodInds2);
        end
    end


    %% Match azimuths and nans

    azRes=round((data1in.azimuth(2)-data1in.azimuth(1))*10)/10;
    if azRes==0.5
        pastDot=data1in.azimuth(1)-floor(data1in.azimuth(1));
        if (pastDot>=0.2 & pastDot<=0.3) | (pastDot>=0.7 & pastDot<=0.8)
            allAz=0.25:azRes:360;
        else
            allAz=0.5:azRes:360;
        end
    else
        allAz=1:360;
    end

    if ~isempty(minMaxAz)
        %         data1in.azimuth(data1in.azimuth<minMaxAz(1) | data1in.azimuth>minMaxAz(2))=nan;
        %         data2in.azimuth(data2in.azimuth<minMaxAz(1) | data2in.azimuth>minMaxAz(2))=nan;
        allAz(allAz<minMaxAz(1))=[];
        allAz(allAz>minMaxAz(2))=[];
    end

    ib1=[];
    ib2=[];
    ibAll=[];
    for kk=1:length(allAz)
        [minDiff1,minInd1]=min(abs(data1in.azimuth-allAz(kk)));
        [minDiff2,minInd2]=min(abs(data2in.azimuth-allAz(kk)));
        if minDiff1<azRes/2 & minDiff2<azRes/2
            ib1=cat(1,ib1,minInd1);
            ib2=cat(1,ib2,minInd2);
            ibAll=cat(1,ibAll,kk);
        end
    end

    data1=[];
    data1.range=data1in.range;
    data2=[];
    data2.range=data2in.range;

    for ii=1:size(inFields,1)
        if ~strcmp(inFields{ii},'range') & ~strcmp(inFields{ii},'time')
            data1.(inFields{ii})=nan(length(allAz),size(data1in.(inFields{ii}),2));
            data1.(inFields{ii})(ibAll,:)=data1in.(inFields{ii})(ib1,:);
            data2.(inFields{ii})=nan(length(allAz),size(data2in.(inFields{ii}),2));
            data2.(inFields{ii})(ibAll,:)=data2in.(inFields{ii})(ib2,:);
        end
    end

    % CMD
    if censorOnCMD
        cmd=zeros(size(data1.DBZ_F));
        if isfield(data1in,'CMD_FLAG')
            data1in.CMD_FLAG=data1in.CMD_FLAG(:,goodInds1);
            cmd(ibAll,:)=data1in.CMD_FLAG(ib1,:);
        elseif isfield(data2in,'CMD_FLAG')
            data2in.CMD_FLAG=data2in.CMD_FLAG(:,goodInds2);
            cmd(ibAll,:)=data2in.CMD_FLAG(ib2,:);
        end
        if sum(cmd,'all')==0
            censorOnCMD=0;
            warning('No CMD flag found.')
        end
    end

    % SNR
    if ~isempty(censorOnSNR)
        snr=zeros(size(data1.DBZ_F));
        if isfield(data1in,'SNR')
            data1in.SNR=data1in.SNR(:,goodInds1);
            snr(ibAll,:)=data1in.SNR(ib1,:);
        elseif isfield(data2in,'SNR')
            data2in.SNR=data2in.SNR(:,goodInds2);
            snr(ibAll,:)=data2in.SNR(ib2,:);
        end
        if isempty(snr)
            censorOnSNR=0;
            disp('No SNR found.')
        end
    end

    % TRIP
    if censorOnTRIP
        trip=zeros(size(data1.DBZ_F));
        if isfield(data1in,'TRIP')
            data1in.TRIP=data1in.TRIP(:,goodInds1);
            trip(ibAll,:)=data1in.TRIP(ib1,:);
        elseif isfield(data2in,'TRIP')
            data2in.TRIP=data2in.TRIP(:,goodInds2);
            trip(ibAll,:)=data2in.TRIP(ib2,:);
        end
        if isempty(trip)
            censorOnTRIP=0;
            disp('No TRIP flag found.')
        end
    end

    for ii=1:size(inFields,1)
        if ~strcmp(inFields{ii},'range') & ~strcmp(inFields{ii},'time')
            % Censor on DBZ
            if censorOnDBZ & size(data1.(inFields{ii}))==size(data1.DBZ_F)
                data1.(inFields{ii})(isnan(data1.DBZ_F))=nan;
                data2.(inFields{ii})(isnan(data2.DBZ_F))=nan;
            end
            % Censor on VEL
            if censorOnVEL & size(data1.(inFields{ii}))==size(data1.DBZ_F)
                data1.(inFields{ii})(data1.VEL_F>-2 & data1.VEL_F<2)=nan;
                data2.(inFields{ii})(isnan(data2.DBZ_F))=nan;
            end
            % Censor on CMD
            if censorOnCMD & size(data1.(inFields{ii}))==size(data1.DBZ_F)
                if ~strcmp(inFields{ii},'CMD_FLAG')
                    data1.(inFields{ii})(cmd<0.5)=nan;
                    data2.(inFields{ii})(cmd<0.5)=nan;
                    data1.(inFields{ii})(isnan(cmd))=nan;
                    data2.(inFields{ii})(isnan(cmd))=nan;
                end
            end
            % Censor on SNR
            if ~isempty(censorOnSNR) & size(data1.(inFields{ii}))==size(data1.DBZ_F)
                data1.(inFields{ii})(snr<censorOnSNR)=nan;
                data2.(inFields{ii})(snr<censorOnSNR)=nan;
            end
            % Censor on TRIP
            if censorOnTRIP & size(data1.(inFields{ii}))==size(data1.DBZ_F)
                data1.(inFields{ii})(trip==1)=nan;
                data2.(inFields{ii})(trip==1)=nan;
            end
            % Match nans
            data1.(inFields{ii})(isnan(data2.(inFields{ii})))=nan;
            data2.(inFields{ii})(isnan(data1.(inFields{ii})))=nan;
        end
    end

    data1.PHIDP_F=wrapTo360(data1.PHIDP_F);
    data2.PHIDP_F=wrapTo360(data2.PHIDP_F);

    %% Plot preparation

    ang_p = deg2rad(90-data1.azimuth);

    angMat=repmat(ang_p,size(data1.range,1),1);

    XX = (data1.range.*cos(angMat));
    YY = (data1.range.*sin(angMat));

    %% Loop through fields
    outstr=inAll{1,3}(aa);
    outstr=outstr{:};

    outDots=strsplit(outstr,'.');
    outParts=strsplit(outDots{1},'VS');

    edges.DBZ_F=[-60:1:60];
    edges.VEL_F=[-50:1:50];
    edges.WIDTH_F=[0:0.2:15];
    edges.PHIDP_F=[0:5:360];
    edges.RHOHV_F=[0.9:0.001:1.1];
    edges.ZDR_F=[-10:0.1:10];

    close all
    f1 = figure('Position',[200 500 1200 750],'DefaultAxesFontSize',12);
    t = tiledlayout(2,3,'TileSpacing','tight','Padding','tight');
    colormap('jet');

    tile=1;

    for jj=1:length(inFields)

        if strcmp(inFields{jj},'azimuth') | strcmp(inFields{jj},'elevation') | ...
                strcmp(inFields{jj},'range') | strcmp(inFields{jj},'time')  | strcmp(inFields{jj},'CMD_FLAG')
            continue
        end

        thisName=inFields{jj};
        xyData=cat(2,data1.(thisName)(:),data2.(thisName)(:));
        xyData(any(isnan(xyData),2),:)=[];
        [N,~,~]=histcounts2(xyData(:,1),xyData(:,2),edges.(thisName),edges.(thisName));
        % N=cat(1,N,N(end,:));
        % N=cat(2,N,N(:,end));

        xyData(any(xyData<edges.(thisName)(1),2),:)=[];
        xyData(any(xyData>edges.(thisName)(end),2),:)=[];


        plotCoords=edges.(thisName)(1:end-1)+(edges.(thisName)(2:end)-edges.(thisName)(1:end-1))/2;
        %plotCoords=cat(2,plotCoords,plotCoords(end)+1);

        s1=nexttile(tile);
        hold on
        if sum(N,'all')>0
            surf(plotCoords,plotCoords,log(N),'edgecolor','none');
            view(2)
            %shading('flat');
            colorbar
            floorAx=prctile(xyData(:,1),1);
            ceilAx=prctile(xyData(:,1),99);
            xlim([floorAx,ceilAx]);
            ylim([floorAx,ceilAx]);
            box on
            xlabel(outParts{2})
            ylabel(outParts{1})
            set(gca,'layer','top')
            xlims=s1.XLim;
            ylims=s1.YLim;
            plot(xlims,ylims,'-k','LineWidth',2);
            s1.SortMethod='childorder';
        end
        title([thisName],'Interpreter','none');
        tile=tile+1;
    end
    set(gcf,'PaperPositionMode','auto')
    print([figdir,outstr,'_2Dhist.png'],'-dpng','-r0');

end