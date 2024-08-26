% Plot hit miss table

clear all
close all

addpath(genpath('~/git/lrose-test/bomb_snowstorm/analysis/'));

figdir='/scr/cirrus1/rsfdata/projects/bomb_snowstorm/figures/paper2024_2/';

indir='/scr/sci/romatsch/forJohn/ams22plots/';

infileList={'35-4svel.txt';
    '37-4svel.txt';
    '39-4svel.txt';
    '41-4svel.txt';
    '43-4svel.txt';
    '45-4svel.txt';
    '47-4svel.txt';
    '49-4svel.txt';
    'svel.dat'};

titles={'(a) Order 35';
    '(b) Order 37';
    '(c) Order 39';
    '(d) Order 41';
    '(e) Order 43';
    '(f) Order 45';
    '(g) Order 47';
    '(h) Order 49';
    '(i) Legacy'};

cm=turbo(32);

figure('Position',[200 500 850 850],'DefaultAxesFontSize',12);
colormap(flipud(cm));
t = tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

xtickLoc=2:2:16;
xtickLab={'1','2','3','4','5','6','7','8'};

ytickLoc=1:5:26;
ytickLab={'0','10','20','30','40','50'};

for kk=1:length(infileList)
    infile=infileList{kk};

indata=table2array(readtable([indir,infile]));

s=nexttile(kk);
hold on

imagesc(flipud(indata));
set(gca,'YDir','normal');
set(gca,'Xtick',xtickLoc);
set(gca,'XtickLabel',xtickLab);
set(gca,'Ytick',ytickLoc);
set(gca,'YtickLabel',ytickLab);
clim([1 5]);

if kk==6
    cb=colorbar;
    cb.Title.String='m s^{-1}';
end

xline=0.5:1:18;
yline=0.5:1:52;

for ii=1:length(xline)
    plot([xline(ii),xline(ii)],[0,27],'-k','linewidth',0.5);
end

for ii=1:length(yline)
    plot([0,17],[yline(ii),yline(ii)],'-k','linewidth',0.5);
end

ax=gca;
ax.SortMethod = 'childorder';

xlim([0.5,16.5]);
ylim([0.5,26.5]);

if kk==7 | kk==8 | kk==9
xlabel('W_1 (m s^{-1})')
end
if kk==1 | kk==4 | kk==7
ylabel('P_1/P_2 (dB)');
end

title(titles{kk})
% ax.Position=[0.123,0.14,0.75,0.84];
end
cb.Layout.Tile = 'east';
%cb.Position=[0.945,0.3,0.02,0.6];

set(gcf,'PaperPositionMode','auto')
print([figdir,'figure4.png'],'-dpng','-r0')
