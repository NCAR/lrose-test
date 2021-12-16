% Compare convStrat output with lightning data

clear all;
close all;

addpath(genpath('~/git/lrose-test/convstrat/dataProcessing/'));

indir=['/scr/cirrus1/rsfdata/projects/convstrat/analysis/gpmPECAN/'];
figdir=['/scr/cirrus1/rsfdata/projects/convstrat/analysis/paperFigs/'];

%% Load data

load([indir,'gpm_20150601_0000_to_20150717_0000.mat']);

%% Calc

% Stratiform
stratPart=partOrigNew;
stratPart(partOrigNew(:,1)~=1,:)=[];

convPart=partOrigNew;
convPart(partOrigNew(:,1)~=2,:)=[];

otherPart=partOrigNew;
otherPart(partOrigNew(:,1)~=3,:)=[];

edges=[13,15,17,20,30,33,35,37,39];
[catCountsStrat]=histcounts(stratPart(:,2),edges);
[catCountsConv]=histcounts(convPart(:,2),edges);
[catCountsOther]=histcounts(otherPart(:,2),edges);

%% Plot
close all

fig=figure('Position',[200 500 900 800],'DefaultAxesFontSize',12);

categories={'Strat low','Strat mid','Strat high','Mixed',...
    'Conv elev','Conv shallow','Conv mid','Conv deep'};

colmapPart=[0,0.1,0.6;
    0.38,0.42,0.96;
    0.65,0.74,0.86;
    0.32,0.78,0.59;
    1,0,1;
    1,1,0;
    0.99,0.77,0.22;
    1 0 0];

close all

wi=5;
hi=9;

fig1=figure('DefaultAxesFontSize',11,'DefaultFigurePaperType','<custom>','units','inch','position',[3,100,wi,hi]);
fig1.PaperPositionMode = 'manual';
fig1.PaperUnits = 'inches';
fig1.Units = 'inches';
fig1.PaperPosition = [0, 0, wi, hi];
fig1.PaperSize = [wi, hi];
fig1.Resize = 'off';
fig1.InvertHardcopy = 'off';

set(fig1,'color','w');

s1=subplot(3,1,1);

b=bar(catCountsStrat./size(stratPart,1).*100,1,'FaceColor','flat');
set(gca,'XTickLabel',categories);
set(gca,'YTick',0:10:100);
for kk = 1:8
    b.CData(kk,:)=colmapPart(kk,:);
end
xtickangle(45);
xlim([0.5,8.5]);
%ylim([0,35]);
set(gca, 'YGrid', 'on');

ylabel('Percent of data points (%)')
title('Stratiform');
box on

s2=subplot(3,1,2);

b=bar(catCountsConv./size(convPart,1).*100,1,'FaceColor','flat');
set(gca,'XTickLabel',categories);
set(gca,'YTick',0:10:100);
for kk = 1:8
    b.CData(kk,:)=colmapPart(kk,:);
end
xtickangle(45);
xlim([0.5,8.5]);
%ylim([0,35]);
set(gca, 'YGrid', 'on');

ylabel('Percent of data points (%)')
title('Convective');
box on

s2=subplot(3,1,3);

b=bar(catCountsOther./size(otherPart,1).*100,1,'FaceColor','flat');
set(gca,'XTickLabel',categories);
set(gca,'YTick',0:10:100);
for kk = 1:8
    b.CData(kk,:)=colmapPart(kk,:);
end
xtickangle(45);
xlim([0.5,8.5]);
%ylim([0,35]);
set(gca, 'YGrid', 'on');

ylabel('Percent of data points (%)')
title('Other');
box on

print([figdir,'gpmVSgpm.png'],'-dpng','-r0')

