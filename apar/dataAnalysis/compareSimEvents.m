% Compare APAR simulated radar data with simulator data

clear all;
close all;

% This command loads the utilities folder
addpath(genpath('~/git/lrose-test/apar/dataAnalysis/utils/'));

%% Specify parameters
% User parameters can be specified in this block. The rest of the script
% should generally be left alone.

event='hurricane';
%event='squall_line';
%event='supercell';

% Display and save plot or just save it
showPlot='off'; % If 'on', plots will be displayed and saved, if 'off' they will only be saved

% Specify minimum range (to ignore erroneous data near the radar)
minRange=0; % Minimum range in km

% Data directories for the input data
indirRad=['/scr/virga1/rsfdata/projects/apar/events/',event,'/sim-ts_thru_rsp_moments/'];
indirSim=['/scr/virga1/rsfdata/projects/apar/events/',event,'/truth_moments/'];

% Ouptut directory for the figures
figdir=['/scr/virga1/rsfdata/projects/apar/events/figures/',event,'/'];

% Specify lower and upper limits of the color scales for each variable
colLims.DBZ=[-10,65];
colLims.VEL=[-30,30];
colLims.WIDTH=[0,5];
colLims.ZDR=[-5,5];

% Specify plot frequency
% When running over a long(ish) time period, it is desirable not to plot every
% single file. The plot frequency can be specified below. E.g. plotFreq=5
% will plot every 5th file. The overall statistics plot will include all
% data not matter what is specified here.
plotFreq=1;

%% List files in input directories
fileListRad=dir([indirRad,'*.nc']);
fileListSim=dir([indirSim,'*.nc']);

%% Sort Sim files into RHIs and PPIs
timeS_PPI=[];
timeS_RHI=[];
indS_PPI=[];
indS_RHI=[];

for bb=1:size(fileListSim,1)
    thisFile=[indirSim,fileListSim(bb).name];
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

%% Loop through radar files
ii=0; % ii keeps track of the files where a match is found
plotInds=1:plotFreq:size(fileListRad,1);

% Initialize comparison variables
dataComp.DBZ=[];
dataComp.VEL=[];
dataComp.WIDTH=[];
dataComp.ZDR=[];

for aa=1:size(fileListRad,1) % Loop through each file

    %% Read radar data
    infileR=[indirRad,fileListRad(aa).name];
    
    disp(['File ',num2str(aa), ' of ',num2str(size(fileListRad,1))]);
    disp(infileR);

    dataR=[];

    dataR.DBZ=[];
    dataR.VEL=[];
    dataR.WIDTH=[];
    dataR.ZDR=[];

    dataFields=fields(dataR);

    dataS=dataR;

    try
        dataR=read_apar(infileR,dataR);
    catch
        continue
    end

    %% Find matching simulator file
    timeR=dataR(1).time(1);

    if strcmp(dataR(1).sweepMode,'sector') % PPI
        [minDiff,sInd]=min(abs(timeS_PPI-timeR));
        if minDiff<seconds(5)
            infileS=[indirSim,fileListSim(indS_PPI(sInd)).name];
        else
            disp('No matching file found.')
            continue
        end
        for kk=1:size(dataR,2)
            sweepM=dataR(kk).sweepMode;
            if ~strcmp(sweepM,'sector')
                dataR(kk)=[];
            end
        end
    elseif strcmp(dataR(1).sweepMode,'rhi') % RHI
        [minDiff,sInd]=min(abs(timeS_RHI-timeR));
        if minDiff<seconds(5)
            infileS=[indirSim,fileListSim(indS_RHI(sInd)).name];
        else
            disp('No matching file found.')
            continue
        end
        for kk=1:size(dataR,2)
            sweepM=dataR(kk).sweepMode;
            if ~strcmp(sweepM,'rhi')
                dataR(kk)=[];
            end
        end
    else
        continue
    end

    % Read simulator data
    dataS=read_apar(infileS,dataS);

    % Check if the size of the two input data sets match
    if ~isequal(size(dataR),size(dataS))
        continue
    end

    %% Loop through scans

    dataVars=fields(dataR);

    for cc=1:size(dataR,2)
        % There are some strange files that have only one ray which we
        % ignore.
        if length(dataR(cc).azimuth)==1
            continue
        end
        if strcmp(dataR(1).sweepMode,'rhi') % For RHIs, check azimuths
            % Check that azimuths match
            if abs(median(dataR(cc).azimuth)-median(dataS(cc).azimuth))>0.1
                disp('Azimuths do not match')
                continue
            end
            angR=round(dataR(cc).elevation,1);
            angS=round(dataS(cc).elevation,1);
        elseif strcmp(dataR(1).sweepMode,'sector') % For PPIs, check elevations
            % Check that elevations match
            if abs(median(dataR(cc).elevation)-median(dataS(cc).elevation))>0.1
                disp('Elevations do not match')
                continue
            end
            angR=round(dataR(cc).azimuth,1);
            angS=round(dataS(cc).azimuth,1);
        end

        % Find matching azimuth or elevation indices
        [bothAngs,ia,ib]=intersect(angR,angS);

        % Trimm both data sources to matching angles
        for dd=1:length(dataVars)
            if ~strcmp(dataVars{dd},'range') & ~strcmp(dataVars{dd},'sweepMode')
                dataRinds(cc).(dataVars{dd})=dataR(cc).(dataVars{dd})(ia,:);
                dataSinds(cc).(dataVars{dd})=dataS(cc).(dataVars{dd})(ib,:);
            end
        end

        dataRinds(cc).range=dataR.range;
        dataSinds(cc).range=dataS.range;

        % Match ranges and missing values
        [bothRange,ia,ib]=intersect(dataRinds(cc).range,dataSinds(cc).range);
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

    % Plot if part of specified plot frequency
    if ismember(ii,plotInds)
        %% Plot preparation
        
        if strcmp(dataR(1).sweepMode,'sector')
            [phi_plt,r2] = meshgrid(deg2rad(dataRinds(1).azimuth),dataRinds(1).range);
            xlimits=[-50,50];
            ylimits=[0,70];
        elseif strcmp(dataR(1).sweepMode,'rhi')
            [phi_plt,r2] = meshgrid(deg2rad(dataRinds(1).elevation),dataRinds(1).range);
            xlimits=[0,70];
            ylimits=[0,20];
        end

        [X,Y]=pol2cart(phi_plt,r2);
        
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
            view(2)

            xlim(xlimits)
            ylim(ylimits)
            clim(colLims.(dataFields{dd}));
            xlabel('Distance (km)');
            ylabel('Distance, (km)');

            colorbar

            title([(dataFields{dd}),' truth ',datestr(fileTime,'yyyy-mm-dd HH:MM:SS')]);

            grid on
            box on

            s3=nexttile(3);

            p=surf(X,Y,dataRinds(1).(dataFields{dd})'-dataSinds(1).(dataFields{dd})', 'EdgeColor','none');
            view(2)

            xlim(xlimits)
            ylim(ylimits)
            clim([-10,10]);
            xlabel('Distance (km)');
            ylabel('Distance, (km)');

            s3.Colormap=redblue;
            colorbar

            title([(dataFields{dd}),' radar-truth ',datestr(fileTime,'yyyy-mm-dd HH:MM:SS')]);

            grid on
            box on

            s4=nexttile(4);

            hold on
            scatter(dataRinds(1).(dataFields{dd})(:),dataSinds(1).(dataFields{dd})(:));
            axis('equal')

            xlabel('Radar');
            ylabel('Truth');

            title([(dataFields{dd}),' radar vs truth ',datestr(fileTime,'yyyy-mm-dd HH:MM:SS')]);

            grid on
            box on

            xlimG=s4.XLim;
            ylimG=s4.YLim;

            minmin=min([xlimG,ylimG]);
            maxmax=max([xlimG,ylimG]);

            plot([minmin,maxmax],[minmin,maxmax],'-','LineWidth',2,'Color',[0.7,0.7,0.7]);

            s4.XLim=xlimG;
            s4.YLim=ylimG;

            print([figdir,'vars/',(dataFields{dd}),'_',datestr(fileTime,'yyyymmdd_HHMMSS'),'.png'],'-dpng','-r0');
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
    ylabel('Truth');

    legend({'Orthogonal fit','1:1'},'Location','northwest');

    grid on
    box on

    title(dataFields{dd});

    ss=ss+1;
end

print([figdir,'scatterStats.png'],'-dpng','-r0');