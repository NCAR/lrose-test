% Read and diplay radar data

clear all;
close all;

addpath(genpath('~/git/lrose-test/apar/dataAnalysis/utils/'));

showPlot='on';

whichPulse='long'; % long or short
whichFile='interleaved'; % interleaved or non_interleaved

if strcmp(whichFile,'interleaved')
    startTime=datetime(2024,5,27,21,31,0);
    endTime=datetime(2024,5,27,21,33,30);
elseif strcmp(whichFile,'non_interleaved')
    startTime=datetime(2024,5,27,21,35,0);
    endTime=datetime(2024,5,27,21,40,30);
end

minRange=6; % Minimum range in km

indirRad=['/scr/virga1/rsfdata/projects/apar/sim/cfradial/moments/',whichFile,'/',whichPulse,'_pulse/'];
indirSim=['/scr/virga1/rsfdata/projects/apar/sim/cfradial/moments/',whichFile,'/aesa_simulator/'];

figdir=['/scr/virga1/rsfdata/projects/apar/sim/cfradial/moments/',whichFile,'/figures/'];

fileListRad=makeFileList(indirRad,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
fileListSim=makeFileList(indirSim,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

% Sort Sim files into RHIs and PPIs
timeS_PPI=[];
timeS_RHI=[];
indS_PPI=[];
indS_RHI=[];

for bb=1:size(fileListSim,1)
    thisFile=fileListSim{bb};
    thisSplit=split(thisFile,["."]);
    if strcmp(thisFile(end-5:end-3),'PPI')
        timeS_PPI=[timeS_PPI;datetime(str2num(thisSplit{2}(1:4)),str2num(thisSplit{2}(5:6)),str2num(thisSplit{2}(7:8)), ...
            str2num(thisSplit{2}(10:11)),str2num(thisSplit{2}(12:13)),str2num(thisSplit{2}(14:15)))];
        indS_PPI=[indS_PPI;bb];
    elseif strcmp(thisFile(end-5:end-3),'RHI')
        timeS_RHI=[timeS_RHI;datetime(str2num(thisSplit{2}(1:4)),str2num(thisSplit{2}(5:6)),str2num(thisSplit{2}(7:8)), ...
            str2num(thisSplit{2}(10:11)),str2num(thisSplit{2}(12:13)),str2num(thisSplit{2}(14:15)))];
        indS_RHI=[indS_RHI;bb];
    end
end

dataComp.DBZ=[];
dataComp.VEL=[];
dataComp.WIDTH=[];
dataComp.ZDR=[];

colLims.DBZ=[-10,65];
colLims.VEL=[-50,50];
colLims.WIDTH=[0,5];
colLims.ZDR=[-5,5];

ii=0;

plotInds=1:5:size(fileListRad,1);

%% Loop through radar files
for aa=1:size(fileListRad,1)

    %% Read radar
    infileR=fileListRad{aa};
    
    disp(['File ',num2str(aa), ' of ',num2str(size(fileListRad,1))]);
    disp(infileR);

    dataR=[];

    dataR.DBZ=[];
    dataR.VEL=[];
    dataR.WIDTH=[];
    dataR.ZDR=[];

    dataFields=fields(dataR);

    dataS=dataR;

    dataR=read_apar(infileR,dataR);

    %% Find matching simulator file
    timeR=dataR(1).time(1);

    if strcmp(dataR(1).sweepMode,'sector')
        [minDiff,sInd]=min(abs(timeS_PPI-timeR));
        if minDiff<seconds(5)
            infileS=fileListSim{indS_PPI(sInd)};
        else
            disp('No matching file found.')
            continue
        end
    % elseif strcmp(dataR(1).sweepMode,'rhi') & strcmp(whichFile,'non_interleaved')
    %     [minDiff,sInd]=min(abs(timeS_RHI-timeR));
    %     if minDiff<seconds(5)
    %         infileS=fileListSim{indS_RHI(sInd)};
    %     else
    %         disp('No matching file found.')
    %         continue
    %     end
    else
        continue
    end

    dataS=read_apar(infileS,dataS);

    %% Loop through scans

    dataVars=fields(dataR);

    for cc=1:size(dataR,2)
        % Round azimuth and elevation
        if length(dataR(cc).azimuth)==1
            continue
        end
        if strcmp(dataR(1).sweepMode,'rhi')
            % Check that azimuths match
            if abs(median(dataR(cc).azimuth)-median(dataS(cc).azimuth))>0.1
                disp('Azimuths do not match')
                continue
            end
            angR=round(dataR(cc).elevation,1);
            angS=round(dataS(cc).elevation,1);
        elseif strcmp(dataR(1).sweepMode,'sector')
            % Check that elevations match
            if abs(median(dataR(cc).elevation)-median(dataS(cc).elevation))>0.1
                disp('Elevations do not match')
                continue
            end
            angR=round(dataR(cc).azimuth,1);
            angS=round(dataS(cc).azimuth,1);
        end

        % Find matching azimuth or elevation inds
        [bothAngs,ia,ib]=intersect(angR,angS);

        for dd=1:length(dataVars)
            if ~strcmp(dataVars{dd},'range') & ~strcmp(dataVars{dd},'sweepMode')
                dataRinds.(dataVars{dd})=dataR(cc).(dataVars{dd})(ia,:);
                dataSinds.(dataVars{dd})=dataS(cc).(dataVars{dd})(ib,:);
            end
        end

        dataRinds.range=dataR.range;
        dataSinds.range=dataS.range;

        % Match range and nans
        [bothRange,ia,ib]=intersect(dataRinds.range,dataSinds.range);
        for dd=1:length(dataFields)
            dataRinds(cc).(dataFields{dd})=dataRinds(cc).(dataFields{dd})(:,ia);
            dataSinds(cc).(dataFields{dd})=dataSinds(cc).(dataFields{dd})(:,ib);

            dataRinds(cc).range=dataRinds(cc).range(ia);
            dataSinds(cc).(dataFields{dd})=dataSinds(cc).(dataFields{dd})(:,ib);
            
            dataRinds(cc).(dataFields{dd})(isnan(dataSinds(cc).(dataFields{dd})))=nan;
            dataSinds(cc).(dataFields{dd})(isnan(dataRinds(cc).(dataFields{dd})))=nan;

            % Trim range
            minRangeInd=max(find(dataRinds(cc).range<minRange));
            dataRinds(cc).(dataFields{dd})(:,1:minRangeInd)=nan;
            dataSinds(cc).(dataFields{dd})(:,1:minRangeInd)=nan;

            vecR=dataRinds(cc).(dataFields{dd})(:);
            vecR(isnan(vecR))=[];
            vecS=dataSinds(cc).(dataFields{dd})(:);
            vecS(isnan(vecS))=[];

            vecRS=cat(2,vecR,vecS);

            dataComp.(dataFields{dd})=cat(1,dataComp.(dataFields{dd}),vecRS);
        end
    end

    ii=ii+1;

    if ismember(ii,plotInds)
        %% Plot preparation
        plotRange=75;
        if strcmp(dataR(1).sweepMode,'sector')
            [phi_plt,r2] = meshgrid(deg2rad(dataRinds(1).azimuth),dataRinds(1).range);
        elseif strcmp(dataR(1).sweepMode,'rhi')
            [phi_plt,r2] = meshgrid(deg2rad(dataRinds(1).elevation),dataRinds(1).range);
        end

        [X,Y] = pol2cart(phi_plt,r2);
        xlimits=[floor(min(X(:))),ceil(max(X(:)))];
        ylimits=[floor(min(Y(:))),ceil(max(Y(:)))];

        for dd=1:length(dataFields)
            if all(isnan(dataRinds(1).(dataFields{dd})(:)))
                continue
            end
            close all

            figure('Position',[200 500 1200 800],'DefaultAxesFontSize',12,'visible',showPlot);
            colormap('jet');

            t = tiledlayout(2,2,'TileSpacing','tight','Padding','tight');

            s1=nexttile(1);

            p=surf(X,Y,dataRinds(1).(dataFields{dd})', 'EdgeColor','none');
            axis('equal')
            view(2)

            xlim(xlimits)
            ylim(ylimits)
            clim(colLims.(dataFields{dd}));
            xlabel('Distance (km)');
            ylabel('Distance, (km)');

            colorbar

            fileTime=mean(dataRinds(1).time);
            title([(dataFields{dd}),' radar ',datestr(fileTime,'yyyy-mm-dd HH:MM:SS')]);

            grid on
            box on

            s2=nexttile(2);

            p=surf(X,Y,dataSinds(1).(dataFields{dd})', 'EdgeColor','none');
            axis('equal')
            view(2)

            xlim(xlimits)
            ylim(ylimits)
            clim(colLims.(dataFields{dd}));
            xlabel('Distance (km)');
            ylabel('Distance, (km)');

            colorbar

            title([(dataFields{dd}),' sim ',datestr(fileTime,'yyyy-mm-dd HH:MM:SS')]);

            grid on
            box on

            s3=nexttile(3);

            p=surf(X,Y,dataRinds(1).(dataFields{dd})'-dataSinds(1).(dataFields{dd})', 'EdgeColor','none');
            axis('equal')
            view(2)

            xlim(xlimits)
            ylim(ylimits)
            clim([-10,10]);
            xlabel('Distance (km)');
            ylabel('Distance, (km)');

            colorbar

            title([(dataFields{dd}),' radar-sim ',datestr(fileTime,'yyyy-mm-dd HH:MM:SS')]);

            grid on
            box on

            s4=nexttile(4);

            hold on
            scatter(dataRinds(1).(dataFields{dd})(:),dataSinds(1).(dataFields{dd})(:));
            axis('equal')

            xlabel('Radar');
            ylabel('Sim');

            title([(dataFields{dd}),' radar vs sim ',datestr(fileTime,'yyyy-mm-dd HH:MM:SS')]);

            grid on
            box on

            xlimG=s4.XLim;
            ylimG=s4.YLim;

            minmin=min([xlimG,ylimG]);
            maxmax=max([xlimG,ylimG]);

            plot([minmin,maxmax],[minmin,maxmax],'-','LineWidth',2,'Color',[0.7,0.7,0.7]);

            s4.XLim=xlimG;
            s4.YLim=ylimG;

            print([figdir,'vars/',(dataFields{dd}),'_',datestr(fileTime,'yyyymmdd_HHMMSS'),'_',whichPulse,'.png'],'-dpng','-r0');
        end
    end
end

disp([num2str(ii),' valid files processed.']);

%% Plot statistics

close all

if isempty(dataComp.VEL)
    figure('Position',[200 500 1200 600],'DefaultAxesFontSize',12,'visible',showPlot);
    colormap('jet');
    t = tiledlayout(1,2,'TileSpacing','tight','Padding','tight');
else
    figure('Position',[200 500 1200 1200],'DefaultAxesFontSize',12,'visible',showPlot);
    colormap('jet');
    t = tiledlayout(2,2,'TileSpacing','tight','Padding','tight');
end

ss=1;
for dd=1:length(dataFields)
    if isempty(dataComp.(dataFields{dd}))
        continue
    end

    % Histcounts
    minField=min(dataComp.(dataFields{dd})(:));
    maxField=max(dataComp.(dataFields{dd})(:));

    gridSp=(maxField-minField)/100;

    fieldEdges.(dataFields{dd})=floor(minField):gridSp:ceil(maxField);
    fieldX.(dataFields{dd})=fieldEdges.(dataFields{dd})(1:end-1)+(fieldEdges.(dataFields{dd})(2)-fieldEdges.(dataFields{dd})(1))/2;

    N.(dataFields{dd})=histcounts2(dataComp.(dataFields{dd})(:,1),dataComp.(dataFields{dd})(:,2),fieldEdges.(dataFields{dd}),fieldEdges.(dataFields{dd}));
    N.(dataFields{dd})(N.(dataFields{dd})==0)=nan;

    % Orthogonal regression
    fitOrth1=gmregress(dataComp.(dataFields{dd})(:,1),dataComp.(dataFields{dd})(:,2));
    fitAll1=[fitOrth1(2) fitOrth1(1)];

    yFit1 = polyval(fitAll1, fieldEdges.(dataFields{dd}));


    % Plot
    s=nexttile(ss);
    hold on
    h=imagesc(fieldX.(dataFields{dd}),fieldX.(dataFields{dd}),N.(dataFields{dd})');
    set(h,'alphadata',~isnan(N.(dataFields{dd})'))
    set(gca,'YDir','normal');

    plot(fieldEdges.(dataFields{dd}), yFit1,'-r','linewidth',2);
    plot([fieldEdges.(dataFields{dd})(1),fieldEdges.(dataFields{dd})(end)], ...
        [fieldEdges.(dataFields{dd})(1),fieldEdges.(dataFields{dd})(end)],'-','LineWidth',2,'Color',[0.7,0.7,0.7]);
    
    xlim([fieldEdges.(dataFields{dd})(1),fieldEdges.(dataFields{dd})(end)]);
    ylim([fieldEdges.(dataFields{dd})(1),fieldEdges.(dataFields{dd})(end)]);

    colormap('turbo');
    colorbar

    xlabel('Radar');
    ylabel('Sim');

    legend({'Orthogonal fit','1:1'},'Location','northwest');

    grid on
    box on

    title(dataFields{dd});

    ss=ss+1;
end

print([figdir,'scatterStats_',whichPulse,'.png'],'-dpng','-r0');