% Plot hit miss table

clear all
close all

indir='/scr/sci/romatsch/data/fromJohn/regModel/widthEst_V3/';
figdir='/scr/cirrus1/rsfdata/projects/bomb_snowstorm/figures/paper2024/';

if ~exist(figdir,'dir')
    mkdir(figdir)
end

infilesAll={'mpow.d','spow.d', ...
    'mvel.d','svel.d', ...
    'mwth.d','swth.d';};

colLims=[-20,1,20; % mpow
    0,4,16; % spow
    -2,2,16; % mvel
    0,4,16; % svel
    -2,2,16; % mwth
    0,4,16]; % swth

titlestr={'(a) Reg. mean power bias (dB)',
    '(b) Reg. power SD (dB)',
    '(c) Reg. mean velocity bias (m s^{-1})',
    '(d) Reg. velocity SD (m s^{-1})',
    '(e) Reg. mean width bias (m s^{-1})',
    '(f) Reg. width SD (m s^{-1})'};

colSD=jet(16);
colSD=colSD(3:end,:);
colSD=cat(1,[1,1,1;0.5,0.5,1],colSD);

colPB=jet(64);
colPB=flipud(colPB(9:end,:));
ls01=linspace(0,1,8);
lsx1=linspace(colPB(end,2),1,8);
lsb=cat(2,ls01',lsx1',ones(8,1));
colPB=cat(1,colPB,lsb);
colPB=cat(1,colPB,[1,1,1;1,0.9,0.9;1,0.8,0.8]);

ls01=linspace(0,1,6);
blues=cat(2,ls01',ls01',ones(6,1));
reds=cat(2,ones(6,1),flipud(ls01'),flipud(ls01'));
colRB=cat(1,blues,reds);
colRB=cat(1,[0,0,0.5;0,0,0.75],colRB,[0.75,0,0;0.5,0,0]);
colRB=flipud(colRB);

f1=figure('Position',[200 0 1150 1300],'DefaultAxesFontSize',12);
t = tiledlayout(3,2,'TileSpacing','tight','Padding','tight');

for ii=1:length(infilesAll)

    infile=infilesAll{ii};

    indata=readData_Torres([indir,infile]);

    velAx=indata.tableOut(:,1,1);
    velAx=[velAx;velAx(end)+(velAx(end)-velAx(end-1))];
    widthAx=indata.width;
    widthAx=[widthAx;widthAx(end)+(widthAx(end)-widthAx(end-1))];


    plotData=squeeze(indata.tableOut(:,2,:));
    plotData=cat(1,plotData,plotData(end,:));
    plotData=cat(2,plotData,plotData(:,end));

    s1=nexttile(ii);
    hold on

    xlabel('Velocity (m s^{-1})');

    title([titlestr{ii}]);

    ylabel('Spectrum width (m s^{-1})');

    if ii==2 | ii==4 | ii==6
        s1.Colormap=colSD;
    elseif ii==1
        s1.Colormap=colPB;
    elseif ii==3 | ii==5
        s1.Colormap=colRB;
    end

    surf(velAx,widthAx,plotData');

    xlim([velAx(1),velAx(end)]);
    ylim([widthAx(1),widthAx(end)]);

    xtickAll=velAx(1:end-1)+(velAx(2:end)-velAx(1:end-1));
    sub.XTick=xtickAll(1:10:end);
    sub.XTickLabel=cellfun(@(x) num2str(x,'%.1f'),{velAx(1:10:end-1)},'un',0);
    sub.XTickLabelRotation=0;

    sub.YTick=widthAx(1:end-1)+(widthAx(2:end)-widthAx(1:end-1))/2;
    sub.YTickLabel=cellfun(@num2str,{widthAx(1:end-1)},'un',0);

    clim(colLims(ii,1:2));
    colorbar

    scatter(4,4,70,'k','*','LineWidth',1.5);

    s1.SortMethod='childorder';
end
set(gcf,'PaperPositionMode','auto')
%print([figdir,'figure6.png'],'-dpng','-r0')
exportgraphics(f1,[figdir,'figure6.png'],'Resolution','300');