% Read and diplay radar data

clear all;
close all;

addpath(genpath('~/git/lrose-test/apar/dataAnalysis/utils/'));

showPlot='on';

startTime=datetime(2024,5,27,21,31,0);
endTime=datetime(2024,5,27,21,33,30);

whichPulse='long'; % long or short pulse

indirRad=['/scr/virga1/rsfdata/projects/apar/sim/cfradial/moments/interleaved/',whichPulse,'_pulse/'];
indirSim='/scr/virga1/rsfdata/projects/apar/sim/cfradial/moments/interleaved/aesa_simulator/';

figdir='/scr/virga1/rsfdata/projects/apar/sim/cfradial/moments/figures/';

fileListRad=makeFileList(indirRad,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);
fileListSim=makeFileList(indirSim,startTime,endTime,'xxxxxx20YYMMDDxhhmmss',1);

% Sort Sim files into RHIs and PPIs
timeS_PPI=[];
timeS_RHI=[];
indS_PPI=[];
indS_RHI=[];

for bb=1:size(fileListSim,1)
    thisFile=fileListSim{bb};
    thisSplit=split(thisFile,["_","."]);
    if strcmp(thisFile(end-5:end-3),'PPI')
        timeS_PPI=[timeS_PPI;datetime(str2num(thisSplit{3}(1:4)),str2num(thisSplit{3}(5:6)),str2num(thisSplit{3}(7:8)), ...
            str2num(thisSplit{4}(1:2)),str2num(thisSplit{4}(3:4)),str2num(thisSplit{4}(5:6)))];
        indS_PPI=[indS_PPI;bb];
    elseif strcmp(thisFile(end-5:end-3),'RHI')
        timeS_RHI=[timeS_RHI;datetime(str2num(thisSplit{3}(1:4)),str2num(thisSplit{3}(5:6)),str2num(thisSplit{3}(7:8)), ...
            str2num(thisSplit{4}(1:2)),str2num(thisSplit{4}(3:4)),str2num(thisSplit{4}(5:6)))];
        indS_RHI=[indS_RHI;bb];
    end
end

dataComp.DBZ=[];
dataComp.VEL=[];
dataComp.WIDTH=[];
dataComp.ZDR=[];

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
    elseif strcmp(dataR(1).sweepMode,'rhi')
        % [minDiff,sInd]=min(abs(timeS_RHI-timeR));
        % if minDiff<seconds(5)
        %     infileS=fileListSim{indS_RHI(sInd)};
        % else
            disp('No matching file found.')
            continue
        % end
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
            dataRinds.(dataFields{dd})=dataRinds(cc).(dataFields{dd})(:,ia);
            dataSinds.(dataFields{dd})=dataSinds(cc).(dataFields{dd})(:,ib);
            
            dataRinds.(dataFields{dd})(isnan(dataSinds.(dataFields{dd})))=nan;
            dataSinds.(dataFields{dd})(isnan(dataRinds.(dataFields{dd})))=nan;

            vecR=dataRinds.(dataFields{dd})(:);
            vecR(isnan(vecR))=[];
            vecS=dataSinds.(dataFields{dd})(:);
            vecS(isnan(vecS))=[];

            vecRS=cat(2,vecR,vecS);

            dataComp.(dataFields{dd})=cat(1,dataComp.(dataFields{dd}),vecRS);
        end
    end


    %% Plot preparation

    %ang_p = deg2rad(90-dataR.azimuth);
    % 
    % angMat=repmat(ang_p,size(dataR.range,1),1);
    % 
    % XX = (dataR.range.*cos(angMat));
    % YY = (dataR.range.*sin(angMat));
    % 
    % %% Z
    % close all
    % 
    % figure('Position',[200 500 2800 1200],'DefaultAxesFontSize',12,'visible',showPlot);
    % 
    % s1=subplot(2,4,1);
    % surf(XX,YY,dataR.DBZ,'edgecolor','none');
    % view(2);
    % caxis([-10 65])
    % title('DBZ (dBZ)')
    % xlabel('km');
    % ylabel('km');
    % s1.Colormap=dbz_default2;
    % cb1=colorbar('XTick',-10:3:65);
    % 
    % grid on
    % box on
    % 
    % %% ZDR
    % 
    % s2=subplot(2,4,2);
    % h=surf(XX,YY,dataR.ZDR_F,'edgecolor','none');
    % view(2);
    % title('ZDR (dB)')
    % xlabel('km');
    % ylabel('km');
    % 
    % s2.Colormap=zdr_default;
    % colLims=[-inf,-20,-2,-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1,1.5,2,2.5,3,4,5,6,8,10,15,20,50,99,inf];
    % applyColorScale(h,dataR.ZDR_F,zdr_default,colLims);
    % 
    % grid on
    % box on
    % 
    % %% VEL
    % 
    % s3=subplot(2,4,3);
    % h3=surf(XX,YY,dataR.VEL_F,'edgecolor','none');
    % view(2);
    % title('VEL (m s^{-1})')
    % xlabel('km');
    % ylabel('km');
    % 
    % grid on
    % box on
    % 
    % colLims=[-inf,-30,-26,-21,-17,-13,-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10,13,17,21,26,30,inf];
    % applyColorScale(h3,dataR.VEL_F,vel_default2,colLims);
    % 
    % 
    % %% WIDTH
    % 
    % s7=subplot(2,4,7);
    % h=surf(XX,YY,dataR.WIDTH_F,'edgecolor','none');
    % view(2);
    % title('WIDTH (m s^{-1})')
    % xlabel('km');
    % ylabel('km');
    % 
    % colLims=[-inf,0,0.5,1,1.5,2,2.5,3,4,5,6,7,8,10,12.5,15,20,25,50,inf];
    % applyColorScale(h,dataR.WIDTH_F,width_default,colLims);
    % 
    % grid on
    % box on
    % 
    % 
    % %% Save first zoom
    % 
    % linkaxes([s1,s2,s3,s4,s5,s6,s7,s8],'xy');
    % outstr=inAll{1,2}(aa);
    % outstr=outstr{:};
    % 
    % xlimits1=[inAll{1,3}(aa),inAll{1,4}(aa)];
    % ylimits1=[inAll{1,5}(aa),inAll{1,6}(aa)];
    % 
    % xlim(xlimits1)
    % ylim(ylimits1)
    % daspect(s1,[1 1 1]);
    % daspect(s2,[1 1 1]);
    % daspect(s3,[1 1 1]);
    % daspect(s4,[1 1 1]);
    % daspect(s5,[1 1 1]);
    % daspect(s6,[1 1 1]);
    % daspect(s7,[1 1 1]);
    % daspect(s8,[1 1 1]);
    % 
    % print([figdir,outstr,'_zoom1.png'],'-dpng','-r0');

end

%% Plot
scatter(dataComp.DBZ(:,1),dataComp.DBZ(:,2));

scatter(dataComp.ZDR(:,1),dataComp.ZDR(:,2));

