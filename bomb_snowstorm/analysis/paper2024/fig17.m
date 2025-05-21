% Read and diplay radar data

clear all;
close all;

addpath(genpath('~/git/lrose-test/bomb_snowstorm/analysis/utils/'));

infile1='/scr/cirrus1/rsfdata/projects/bomb_snowstorm/tables/SPOL20190313_220622_INDX_CMD_RHV_GAUSS_WN128_V3.txt';

infile2='/scr/cirrus1/rsfdata/projects/bomb_snowstorm/tables/SPOL20190313_220622_INDX_CMD_RHV_GAUSS_REG_V3.txt';

figdir='/scr/cirrus1/rsfdata/projects/bomb_snowstorm/figures/paper2024/';

xlimits1=[-20,100];
ylimits1=[-20,100];

kernel=[9,5]; % Az and range of std kernel. Default: [9,5]

censorOnCMD=1;
%%%%%%%%%%%%%%
censorOnSNR=[]; % Set to empty if not used !!!!!!! Only use areas with SNR above XX dB
%%%%%%%%%%%%%%
halfNyquist=0; % In some files the nyquist needs to be divided by 2

%% Read data


data1in=readDataTables(infile1,' ');
data1in.azimuth=round(data1in.azimuth);
if isfield(data1in,'TRIP')
    data1in.SNR=data1in.TRIP;
    data1in=rmfield(data1in,'TRIP');
end

data2in=readDataTables(infile2,' ');
%data2in.azimuth=round(data2in.azimuth);
if isfield(data2in,'TRIP')
    data2in.SNR=data2in.TRIP;
    data2in=rmfield(data2in,'TRIP');
end

nyquist=26.675;


%% Cut range
inFields1=fields(data1in);
inFields2=fields(data2in);
inFields=intersect(inFields1,inFields2);

minMaxRange=[min(data1in.range),max(data1in.range)];

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
    if isempty(cmd)
        censorOnCMD=0;
        disp('No CMD flag found.')
    end
end


for ii=1:size(inFields,1)
    if ~strcmp(inFields{ii},'range') & ~strcmp(inFields{ii},'time')
        % Censor on CMD
        if censorOnCMD & size(data1.(inFields{ii}))==size(data1.DBZ_F)
            data1.(inFields{ii})(cmd==0)=nan;
            data2.(inFields{ii})(cmd==0)=nan;
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
for jj=1:length(inFields)

    if strcmp(inFields{jj},'azimuth') | strcmp(inFields{jj},'elevation') | ...
            strcmp(inFields{jj},'range') | strcmp(inFields{jj},'time')
        continue
    end


    %% Standard deviations

    if strcmp(inFields{jj},'VEL_F')
        [stdVar1,~]=fast_nd_std(data1.(inFields{jj}),kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',nyquist);
        [stdVar2,~]=fast_nd_std(data2.(inFields{jj}),kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',nyquist);
    elseif strcmp(inFields{jj},'PHIDP_F')
        [stdVar1,~]=fast_nd_std(data1.(inFields{jj}),kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',180);
        [stdVar2,~]=fast_nd_std(data2.(inFields{jj}),kernel,'mode','partial','nan_std',1,'circ_std',1,'nyq',180);
    else
        [stdVar1,~]=fast_nd_std(data1.(inFields{jj}),kernel,'mode','partial','nan_std',1);
        [stdVar2,~]=fast_nd_std(data2.(inFields{jj}),kernel,'mode','partial','nan_std',1);
    end

    stdVar1(isnan(data1.(inFields{jj})))=nan;
    stdVar2(isnan(data2.(inFields{jj})))=nan;

    stdVar1(stdVar1==Inf)=nan;
    stdVar2(stdVar2==Inf)=nan;

    diffField.(inFields{jj})=stdVar2-stdVar1;

end

%% Plot
close all

f1=figure('Position',[200 500 1000 585],'DefaultAxesFontSize',12);
colormap('jet');
t = tiledlayout(2,3,'TileSpacing','compact','Padding','tight');

s1=nexttile(1);
hold on
spacing=0.1;
edges=-1.5:spacing:1.5;
hc=histcounts(diffField.DBZ_F(:),edges);
bar(edges(1:end-1)+spacing/2,hc/sum(hc)*100,1)
xlim([-1.5,1.5]);

xlabel('St. dev. Reg. - st. dev. WN (dB)');
ylabel('Percent of data points (%)');

xticks(-2:0.5:2);
s1.XTickLabelRotation=0;

ylims=s1.YLim;
plot([0,0],ylims,'-r','LineWidth',2);

s1.SortMethod='childorder';
title('(a) Reflectivity (dB)')

grid on
box on

s2=nexttile(2);
hold on
spacing=0.05;
edges=-0.8:spacing:0.8;
hc=histcounts(diffField.VEL_F(:),edges);
bar(edges(1:end-1)+spacing/2,hc/sum(hc)*100,1)
xlim([-0.8,0.8]);

xlabel('St. dev. Reg. - st. dev. WN (dB)');

xticks(-2:0.25:2);
xtickangle(0);

ylims=s2.YLim;
plot([0,0],ylims,'-r','LineWidth',2);

s2.SortMethod='childorder';
title('(b) Velocity (m s^{-1})')

grid on
box on

s3=nexttile(3);
hold on
spacing=0.2;
edges=-3:spacing:3;
hc=histcounts(diffField.PHIDP_F(:),edges);
bar(edges(1:end-1)+spacing/2,hc/sum(hc)*100,1)
xlim([-3,3]);

xlabel(['St. dev. Reg. - st. dev. WN (',char(176),')']);

xticks(-3:3);
xtickangle(0);

ylims=s3.YLim;
plot([0,0],ylims,'-r','LineWidth',2);

s3.SortMethod='childorder';
title(['(c) \phi_{DP} (',char(176),')'])

grid on
box on

s4=nexttile(4);
hold on
spacing=0.025;
edges=-0.4:spacing:0.4;
hc=histcounts(diffField.ZDR_F(:),edges);
bar(edges(1:end-1)+spacing/2,hc/sum(hc)*100,1)
xlim([-0.4,0.4]);

xlabel('St. dev. Reg. - st. dev. WN (dB)');
ylabel('Percent of data points (%)');

xticks(-2:0.1:2);
xtickangle(0);

ylims=s4.YLim;
plot([0,0],ylims,'-r','LineWidth',2);

s4.SortMethod='childorder';
title('(d) Z_{DR} (dB)')

grid on
box on

s5=nexttile(5);
hold on
spacing=0.025;
edges=-0.4:spacing:0.4;
hc=histcounts(diffField.WIDTH_F(:),edges);
bar(edges(1:end-1)+spacing/2,hc/sum(hc)*100,1)
xlim([-0.4,0.4]);

xlabel('St. dev. Reg. - st. dev. WN (m s^{-1})');

xticks(-2:0.1:2);
xtickangle(0);

ylims=s5.YLim;
plot([0,0],ylims,'-r','LineWidth',2);

s5.SortMethod='childorder';
title('(e) Spectrum width (m s^{-1})')

grid on
box on

s6=nexttile(6);
hold on
spacing=0.001;
edges=-0.015:spacing:0.015;
hc=histcounts(diffField.RHOHV_NNC_F(:),edges);
bar(edges(1:end-1)+spacing/2,hc/sum(hc)*100,1)
xlim([-0.015,0.015]);

xlabel('St. dev. Reg. - st. dev. WN');

xticks(-1:0.005:1);
xtickangle(0);

ylims=s6.YLim;
plot([0,0],ylims,'-r','LineWidth',2);

s6.SortMethod='childorder';
title('(f) \rho_{HV}')

grid on
box on

%print([figdir,'figure20.png'],'-dpng','-r0');
exportgraphics(f1,[figdir,'figure17.png'],'Resolution','300');