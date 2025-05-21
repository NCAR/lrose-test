% Read and diplay radar data

clear all;
close all;

addpath(genpath('~/git/lrose-test/bomb_snowstorm/analysis/utils/'));

figdir='/scr/cirrus1/rsfdata/projects/bomb_snowstorm/figures/paper2024/';

xlimits1=[-120,60];
ylimits1=[-110,70];

kernel=[9,5]; % Az and range of std kernel. Default: [9,5]

censorOnCMD=0;
%%%%%%%%%%%%%%
censorOnSNR=[]; % Set to empty if not used !!!!!!! Only use areas with SNR above XX dB
%%%%%%%%%%%%%%
halfNyquist=0; % In some files the nyquist needs to be divided by 2

%% Read data

infile1='/scr/sleet1/rsfdata/projects/eolbase/cfradial/kftg/moments/20220329/cfrad.20220329_221646.829_to_20220329_222242.984_KFTG_SUR.nc';

data1in=[];

data1in.DBZ=[];
data1in.VEL=[];
data1in.WIDTH=[];
data1in.ZDR=[];
data1in.PHIDP=[];
data1in.RHOHV=[];

data1in=read_spol(infile1,data1in);
nyquist=ncread(infile1,'nyquist_velocity');

data1in=data1in(1);

data1in.DBZ_F=data1in.DBZ;
data1in.VEL_F=data1in.VEL;
data1in.WIDTH_F=data1in.WIDTH;
data1in.ZDR_F=data1in.ZDR;
data1in.PHIDP_F=data1in.PHIDP;
data1in.RHOHV_NNC_F=data1in.RHOHV;

data1in.azimuth=round(data1in.azimuth);
if isfield(data1in,'TRIP')
    data1in.SNR=data1in.TRIP;
    data1in=rmfield(data1in,'TRIP');
end

infile2='/scr/cirrus1/rsfdata/projects/nexrad/tables/KFTG_LPRT_SR_20220329_221645_1.91_253.46_14ptFlt28pt-V4.txt';

data2in=readDataTables(infile2,' ');
%data2in.azimuth=round(data2in.azimuth);
if isfield(data2in,'TRIP')
    data2in.SNR=data2in.TRIP;
    data2in=rmfield(data2in,'TRIP');
end

%% Cut range
inFields1=fields(data1in);
inFields2=fields(data2in);
inFields=intersect(inFields1,inFields2);

minMaxRangeOrig=[];
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

minMaxAz=[];
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
if censorOnCMD
    cmd=zeros(size(data2.DBZ_F));
    if isfield(data1in,'CMD_FLAG')
        data1in.CMD_FLAG=data1in.CMD_FLAG(:,goodInds1);
        cmd(ibAll,:)=data1in.CMD_FLAG(ib1,:);
    elseif isfield(data2in,'CMD_FLAG')
        data2in.CMD_FLAG=data2in.CMD_FLAG(:,goodInds2);
        cmd(ibAll,:)=data2in.CMD_FLAG(ib2,:);
    end
    if isempty(cmd)
        censorOnCMD=0;
        disp('No CMD flag found.')
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

for ii=1:size(inFields,1)
    if ~strcmp(inFields{ii},'range') & ~strcmp(inFields{ii},'time')
        % Censor on CMD
        if censorOnCMD & size(data1.(inFields{ii}))==size(data1.DBZ_F)
            data1.(inFields{ii})(cmd==0)=nan;
            data2.(inFields{ii})(cmd==0)=nan;
        end
        % Censor on SNR
        if ~isempty(censorOnSNR) & size(data1.(inFields{ii}))==size(data1.DBZ_F)
            data1.(inFields{ii})(snr<censorOnSNR)=nan;
            data2.(inFields{ii})(snr<censorOnSNR)=nan;
        end
        % Match nans
        data1.(inFields{ii})(isnan(data2.(inFields{ii})))=nan;
        data2.(inFields{ii})(isnan(data1.(inFields{ii})))=nan;
    end
end

%% Plot preparation

ang_p = deg2rad(90-data1.azimuth);

angMat=repmat(ang_p,size(data1.range,1),1);

XX = (data1.range.*cos(angMat));
YY = (data1.range.*sin(angMat));

%% Loop through fields

jj=1;
inFields{jj}='ZDR_F';

%% Standard deviations

if strcmp(inFields{jj},'VEL_F')
    [stdVar1_1,~]=fast_nd_std(data1.(inFields{jj}),kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',nyquist);
    [stdVar2_1,~]=fast_nd_std(data2.(inFields{jj}),kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',nyquist);
elseif strcmp(inFields{jj},'PHIDP_F')
    [stdVar1_1,~]=fast_nd_std(data1.(inFields{jj}),kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',180);
    [stdVar2_1,~]=fast_nd_std(data2.(inFields{jj}),kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',180);
else
    [stdVar1_1,~]=fast_nd_std(data1.(inFields{jj}),kernel,'mode','partial','nan_std',1);
    [stdVar2_1,~]=fast_nd_std(data2.(inFields{jj}),kernel,'mode','partial','nan_std',1);
end

stdVar1_1(isnan(data1.(inFields{jj})))=nan;
stdVar2_1(isnan(data2.(inFields{jj})))=nan;

stdVar1_1(stdVar1_1==Inf)=nan;
stdVar2_1(stdVar2_1==Inf)=nan;

%% Line 2
censorOnCMD=0;
%%%%%%%%%%%%%%
censorOnSNR=[]; % Set to empty if not used !!!!!!! Only use areas with SNR above XX dB
%%%%%%%%%%%%%%
halfNyquist=0; % In some files the nyquist needs to be divided by 2

%% Read data

infile1='/scr/sleet1/rsfdata/projects/eolbase/cfradial/kftg/moments/20220329/cfrad.20220329_221646.829_to_20220329_222242.984_KFTG_SUR.nc';

data1in=[];

data1in.DBZ=[];
data1in.VEL=[];
data1in.WIDTH=[];
data1in.ZDR=[];
data1in.PHIDP=[];
data1in.RHOHV=[];

data1in=read_spol(infile1,data1in);
nyquist=ncread(infile1,'nyquist_velocity');

data1in=data1in(1);

data1in.DBZ_F=data1in.DBZ;
data1in.VEL_F=data1in.VEL;
data1in.WIDTH_F=data1in.WIDTH;
data1in.ZDR_F=data1in.ZDR;
data1in.PHIDP_F=data1in.PHIDP;
data1in.RHOHV_NNC_F=data1in.RHOHV;

data1in.azimuth=round(data1in.azimuth);
if isfield(data1in,'TRIP')
    data1in.SNR=data1in.TRIP;
    data1in=rmfield(data1in,'TRIP');
end

infile2='/scr/cirrus1/rsfdata/projects/nexrad/tables/KFTG_LPRT_SR_20220329_221645_1.91_253.46_22ptFlt28pt-V4.txt';

data2in=readDataTables(infile2,' ');
%data2in.azimuth=round(data2in.azimuth);
if isfield(data2in,'TRIP')
    data2in.SNR=data2in.TRIP;
    data2in=rmfield(data2in,'TRIP');
end

%% Cut range
inFields1=fields(data1in);
inFields2=fields(data2in);
inFields=intersect(inFields1,inFields2);

minMaxRangeOrig=[];
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

minMaxAz=[];
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
if censorOnCMD
    cmd=zeros(size(data2.DBZ_F));
    if isfield(data1in,'CMD_FLAG')
        data1in.CMD_FLAG=data1in.CMD_FLAG(:,goodInds1);
        cmd(ibAll,:)=data1in.CMD_FLAG(ib1,:);
    elseif isfield(data2in,'CMD_FLAG')
        data2in.CMD_FLAG=data2in.CMD_FLAG(:,goodInds2);
        cmd(ibAll,:)=data2in.CMD_FLAG(ib2,:);
    end
    if isempty(cmd)
        censorOnCMD=0;
        disp('No CMD flag found.')
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

for ii=1:size(inFields,1)
    if ~strcmp(inFields{ii},'range') & ~strcmp(inFields{ii},'time')
        % Censor on CMD
        if censorOnCMD & size(data1.(inFields{ii}))==size(data1.DBZ_F)
            data1.(inFields{ii})(cmd==0)=nan;
            data2.(inFields{ii})(cmd==0)=nan;
        end
        % Censor on SNR
        if ~isempty(censorOnSNR) & size(data1.(inFields{ii}))==size(data1.DBZ_F)
            data1.(inFields{ii})(snr<censorOnSNR)=nan;
            data2.(inFields{ii})(snr<censorOnSNR)=nan;
        end
        % Match nans
        data1.(inFields{ii})(isnan(data2.(inFields{ii})))=nan;
        data2.(inFields{ii})(isnan(data1.(inFields{ii})))=nan;
    end
end

%% Loop through fields

jj=1;
inFields{jj}='ZDR_F';

%% Standard deviations

if strcmp(inFields{jj},'VEL_F')
    [stdVar1_2,~]=fast_nd_std(data1.(inFields{jj}),kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',nyquist);
    [stdVar2_2,~]=fast_nd_std(data2.(inFields{jj}),kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',nyquist);
elseif strcmp(inFields{jj},'PHIDP_F')
    [stdVar1_2,~]=fast_nd_std(data1.(inFields{jj}),kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',180);
    [stdVar2_2,~]=fast_nd_std(data2.(inFields{jj}),kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',180);
else
    [stdVar1_2,~]=fast_nd_std(data1.(inFields{jj}),kernel,'mode','partial','nan_std',1);
    [stdVar2_2,~]=fast_nd_std(data2.(inFields{jj}),kernel,'mode','partial','nan_std',1);
end

stdVar1_2(isnan(data1.(inFields{jj})))=nan;
stdVar2_2(isnan(data2.(inFields{jj})))=nan;

stdVar1_2(stdVar1_2==Inf)=nan;
stdVar2_2(stdVar2_2==Inf)=nan;

%% Line 2

censorOnCMD=1;
%%%%%%%%%%%%%%
censorOnSNR=[]; % Set to empty if not used !!!!!!! Only use areas with SNR above XX dB
%%%%%%%%%%%%%%
halfNyquist=0; % In some files the nyquist needs to be divided by 2

%% Read data

infile1='/scr/sleet1/rsfdata/projects/eolbase/cfradial/kftg/moments/20220329/cfrad.20220329_221646.829_to_20220329_222242.984_KFTG_SUR.nc';

data1in=[];

data1in.DBZ=[];
data1in.VEL=[];
data1in.WIDTH=[];
data1in.ZDR=[];
data1in.PHIDP=[];
data1in.RHOHV=[];

data1in=read_spol(infile1,data1in);
nyquist=ncread(infile1,'nyquist_velocity');

data1in=data1in(1);

data1in.DBZ_F=data1in.DBZ;
data1in.VEL_F=data1in.VEL;
data1in.WIDTH_F=data1in.WIDTH;
data1in.ZDR_F=data1in.ZDR;
data1in.PHIDP_F=data1in.PHIDP;
data1in.RHOHV_NNC_F=data1in.RHOHV;

data1in.azimuth=round(data1in.azimuth);
if isfield(data1in,'TRIP')
    data1in.SNR=data1in.TRIP;
    data1in=rmfield(data1in,'TRIP');
end

infile2='/scr/cirrus1/rsfdata/projects/nexrad/tables/KFTG_LPRT_SR_20220329_221645_1.91_253.46_22ptFlt28pt-V4.txt';

data2in=readDataTables(infile2,' ');
%data2in.azimuth=round(data2in.azimuth);
if isfield(data2in,'TRIP')
    data2in.SNR=data2in.TRIP;
    data2in=rmfield(data2in,'TRIP');
end

%% Cut range
inFields1=fields(data1in);
inFields2=fields(data2in);
inFields=intersect(inFields1,inFields2);

minMaxRangeOrig=[];
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

minMaxAz=[];
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
if censorOnCMD
    cmd=zeros(size(data2.DBZ_F));
    if isfield(data1in,'CMD_FLAG')
        data1in.CMD_FLAG=data1in.CMD_FLAG(:,goodInds1);
        cmd(ibAll,:)=data1in.CMD_FLAG(ib1,:);
    elseif isfield(data2in,'CMD_FLAG')
        data2in.CMD_FLAG=data2in.CMD_FLAG(:,goodInds2);
        cmd(ibAll,:)=data2in.CMD_FLAG(ib2,:);
    end
    if isempty(cmd)
        censorOnCMD=0;
        disp('No CMD flag found.')
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

for ii=1:size(inFields,1)
    if ~strcmp(inFields{ii},'range') & ~strcmp(inFields{ii},'time')
        % Censor on CMD
        if censorOnCMD & size(data1.(inFields{ii}))==size(data1.DBZ_F)
            data1.(inFields{ii})(cmd==0)=nan;
            data2.(inFields{ii})(cmd==0)=nan;
        end
        % Censor on SNR
        if ~isempty(censorOnSNR) & size(data1.(inFields{ii}))==size(data1.DBZ_F)
            data1.(inFields{ii})(snr<censorOnSNR)=nan;
            data2.(inFields{ii})(snr<censorOnSNR)=nan;
        end
        % Match nans
        data1.(inFields{ii})(isnan(data2.(inFields{ii})))=nan;
        data2.(inFields{ii})(isnan(data1.(inFields{ii})))=nan;
    end
end

%% Plot preparation

ang_p = deg2rad(90-data1.azimuth);

angMat=repmat(ang_p,size(data1.range,1),1);

XX = (data1.range.*cos(angMat));
YY = (data1.range.*sin(angMat));

%% Loop through fields

jj=1;
inFields{jj}='ZDR_F';

%% Standard deviations

if strcmp(inFields{jj},'VEL_F')
    [stdVar1_3,~]=fast_nd_std(data1.(inFields{jj}),kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',nyquist);
    [stdVar2_3,~]=fast_nd_std(data2.(inFields{jj}),kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',nyquist);
elseif strcmp(inFields{jj},'PHIDP_F')
    [stdVar1_3,~]=fast_nd_std(data1.(inFields{jj}),kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',180);
    [stdVar2_3,~]=fast_nd_std(data2.(inFields{jj}),kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',180);
else
    [stdVar1_3,~]=fast_nd_std(data1.(inFields{jj}),kernel,'mode','partial','nan_std',1);
    [stdVar2_3,~]=fast_nd_std(data2.(inFields{jj}),kernel,'mode','partial','nan_std',1);
end

stdVar1_3(isnan(data1.(inFields{jj})))=nan;
stdVar2_3(isnan(data2.(inFields{jj})))=nan;

stdVar1_3(stdVar1_3==Inf)=nan;
stdVar2_3(stdVar2_3==Inf)=nan;


%% Plot
close all

f1=figure('Position',[200 500 900 1200],'DefaultAxesFontSize',12);
colormap('jet');
t = tiledlayout(3,2,'TileSpacing','tight','Padding','tight');

s1=nexttile(1);
diffField=stdVar2_1-stdVar1_1;
surf(XX,YY,diffField,'edgecolor','none');
view(2);
clim([-5,5]);
s1.Colormap=velCols;
colorbar;
title('(a) 14 pt, filtered and unfiltered');
xlabel('km');
ylabel('km');

grid on
box on

xlim(xlimits1)
ylim(ylimits1)

s2=nexttile(2);

hold on
spacing=0.05;
edges=-1:spacing:1;
hc=histcounts(diffField(:),edges);
bar(edges(1:end-1)+spacing/2,hc/sum(hc)*100,1)
xlim([-1,1]);

xlabel('St. dev. Reg. - st. dev. Level2 (dB)');
ylabel('Percent of data points (%)');

ylims=s2.YLim;
plot([0,0],ylims,'-r','LineWidth',2);

s2.SortMethod='childorder';
title('(b) 14 pt, filtered and unfiltered')

grid on
box on

s3=nexttile(3);
diffField=stdVar2_2-stdVar1_2;
surf(XX,YY,diffField,'edgecolor','none');
view(2);
clim([-5,5]);
s3.Colormap=velCols;
colorbar;
title('(c) 22 pt, filtered and unfiltered');
xlabel('km');
ylabel('km');

grid on
box on

xlim(xlimits1)
ylim(ylimits1)

s4=nexttile(4);

hold on
spacing=0.05;
edges=-1:spacing:1;
hc=histcounts(diffField(:),edges);
bar(edges(1:end-1)+spacing/2,hc/sum(hc)*100,1)
xlim([-1,1]);

xlabel('St. dev. Reg. - st. dev. Level2 (dB)');
ylabel('Percent of data points (%)');

ylims=s4.YLim;
plot([0,0],ylims,'-r','LineWidth',2);

s4.SortMethod='childorder';
title('(d) 22 pt, filtered and unfiltered')

grid on
box on

s5=nexttile(5);
diffField=stdVar2_3-stdVar1_3;
surf(XX,YY,diffField,'edgecolor','none');
view(2);
clim([-5,5]);
s5.Colormap=velCols;
colorbar;
title('(e) 22 pt, filtered');
xlabel('km');
ylabel('km');

grid on
box on

xlim(xlimits1)
ylim(ylimits1)

s6=nexttile(6);

hold on
spacing=0.05;
edges=-1:spacing:1;
hc=histcounts(diffField(:),edges);
bar(edges(1:end-1)+spacing/2,hc/sum(hc)*100,1)
xlim([-1,1]);

xlabel('St. dev. Reg. - st. dev. Level2 (dB)');
ylabel('Percent of data points (%)');

ylims=s6.YLim;
plot([0,0],ylims,'-r','LineWidth',2);

s6.SortMethod='childorder';
title('(f) 22 pt, filtered')

grid on
box on

daspect(s1,[1 1 1]);
daspect(s3,[1 1 1]);
daspect(s5,[1 1 1]);

%print([figdir,'figure26.png'],'-dpng','-r0');
exportgraphics(f1,[figdir,'figure23.png'],'Resolution','300');