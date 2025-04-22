% Plot hit miss table

clear all
close all

addpath(genpath('~/git/lrose-test/bomb_snowstorm/analysis/'));

figdir='/scr/cirrus1/rsfdata/projects/bomb_snowstorm/figures/paper2024_2/';

indir='/scr/sci/romatsch/forJohn/ams22plots/';

infileList={'37-4svel-10SNR.txt';
    'svel-40ptNW-W4-SNR10.txt';
    'svel-48ptNW-W4-SNR10.txt';
    'svel-48ptNW-W2-SNR10.txt'};

cm=turbo(32);

figure('Position',[200 500 600 570],'DefaultAxesFontSize',12);
colormap(flipud(cm));
t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

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

    if kk==2
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

    if kk==3 | kk==4
        xlabel('W_1 (m s^{-1})')
    end

    if kk==3 | kk==1
    ylabel('P_1/P_2 (dB)');
    end


    title(infile(1:end-4))
    % ax.Position=[0.123,0.14,0.75,0.84];
end
cb.Layout.Tile = 'east';
%cb.Position=[0.945,0.3,0.02,0.6];

set(gcf,'PaperPositionMode','auto')
print([figdir,'figureA2.png'],'-dpng','-r0')
