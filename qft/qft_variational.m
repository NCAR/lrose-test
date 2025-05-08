clear all
close all

tic

%figdir='/scr/virga1/rsfdata/projects/qft/d5_deltaVaries_lambda1/';
figdir='/scr/virga1/rsfdata/projects/qft/dVaries_delta1overD_lambda1/';

% Initialize
d=8;
nMax=10;

w=1:0.25:10;
lambda=1;
delta=1/d;
%delta=0.05:0.05:1;
%m2=-(fliplr(0:2:20));
m2=-(fliplr(12:1:16));

%cutoffAll=100:100:1000;
%cutoffAll=[100,500,1000];
cutoffAll=[100,1500];

if ~exist(figdir,'dir')
    mkdir(figdir)
end

%% Create pregens

% Create combinations
vecIn=(0:nMax)';
matIn=repmat(vecIn,1,d);
cellIn=num2cell(matIn,1);
sub=cell(1,numel(cellIn));
[sub{:}]=ndgrid(cellIn{:});
sub=cellfun(@(x)x(:),sub,'UniformOutput', false);
allPerms = cell2mat(sub);
C=fliplr(allPerms);

% sort and unique
uSorted=sort(C,2);
uC=unique(uSorted,'rows');

sumC=sum(uC,2);
uSumC=cat(2,uC,sumC);

sC=sortrows(uSumC,size(uSumC,2));

% Split into odd and even
evensS=[];
oddsS=[];

for ii=1:size(sC,1)
    if mod(sC(ii,end),2)==0
        evensS=cat(1,evensS,sC(ii,:));
    else
        oddsS=cat(1,oddsS,sC(ii,:));
    end
end

%% Loop over delta

for ll=1:length(delta)

    close all
    f1 = figure('Position',[200 500 600 600],'DefaultAxesFontSize',12);

    t = tiledlayout(1,1,'TileSpacing','tight','Padding','tight');
    s1=nexttile(1);
    colmap=turbo(length(cutoffAll));

    hold on

    %% Loop over cutoff

    bcE=nan(length(cutoffAll),length(m2));
    bcO=nan(length(cutoffAll),length(m2));

    for kk=1:length(cutoffAll)
        cutoff=cutoffAll(kk);

        disp(['Cutoff: ',num2str(cutoff)]);

        minAllE=nan(1,length(m2));
        minAllEw=nan(1,length(m2));
        vecAllE=nan(cutoff,length(m2));
        x2exE=nan(1,length(m2));
        x4exE=nan(1,length(m2));

        minAllO=nan(1,length(m2));
        minAllOw=nan(1,length(m2));
        vecAllO=nan(cutoff,length(m2));
        x2exO=nan(1,length(m2));
        x4exO=nan(1,length(m2));

        disp('Even matrix');
        [matOutE1,matOutE2,matOutE3]=calcMatrix(evensS,d,w',delta(ll),cutoff);
        disp('Odd matrix');
        [matOutO1,matOutO2,matOutO3]=calcMatrix(oddsS,d,w',delta(ll),cutoff);

        for jj=1:length(m2)

            disp(['Loop over m: ',num2str(jj),' of ',num2str(length(m2))]);

            % Evens

            matOutE=matOutE1+matOutE2.*(m2(jj)+2./delta(ll).^2)./2+matOutE3.*lambda./delta(ll);

            % Eigen values
            minE=nan(1,length(w));
            vE=nan(cutoff,length(w));

            for ii=1:length(w)
                [V,D]=eig(matOutE(:,:,ii));
                ee=(D(sub2ind(size(D),1:size(D,1),1:size(D,2))))';
                [minE(ii),minEind]=min(ee);
                vE(:,ii)=V(:,minEind);
            end

            [minAllE(jj),minAllEind]=min(minE);
            minAllEw(jj)=w(minAllEind);
            vecAllE(:,jj)=vE(:,minAllEind);

            matThisE2=matOutE2(:,:,minAllEind);
            x2exE(jj)=vecAllE(:,jj)'*matThisE2*vecAllE(:,jj);
            matThisE3=matOutE3(:,:,minAllEind);
            x4exE(jj)=vecAllE(:,jj)'*matThisE3*vecAllE(:,jj).*d;

            % Odds

            matOutO=matOutO1+matOutO2.*(m2(jj)+2./delta(ll).^2)./2+matOutO3.*lambda./delta(ll);

            % Eigen values
            minO=nan(1,length(w));
            vO=nan(cutoff,length(w));

            for ii=1:length(w)
                [V,D]=eig(matOutO(:,:,ii));
                eo=(D(sub2ind(size(D),1:size(D,1),1:size(D,2))))';
                [minO(ii),minOind]=min(eo);
                vO(:,ii)=V(:,minOind);
            end

            [minAllO(jj),minAllOind]=min(minO);
            minAllOw(jj)=w(minAllOind);
            vecAllO(:,jj)=vO(:,minAllOind);

            matThisO2=matOutO2(:,:,minAllOind);
            x2exO(jj)=vecAllO(:,jj)'*matThisO2*vecAllO(:,jj);
            matThisO3=matOutO3(:,:,minAllOind);
            x4exO(jj)=vecAllO(:,jj)'*matThisO3*vecAllO(:,jj).*d;

        end

        bcE(kk,:)=1-x4exE./x2exE.^2./(3.*delta(ll).*d);
        bcO(kk,:)=1-x4exO./x2exO.^2./(3.*delta(ll).*d);

        plot(m2,bcE(kk,:),'-','LineWidth',2,'Color',colmap(kk,:));
        plot(m2,bcO(kk,:),'--','LineWidth',2,'Color',colmap(kk,:));
        drawnow

        outTable=table(m2',minAllEw',minAllE',x2exE',x4exE',bcE',minAllOw',minAllO',x2exO',x4exO',bcO', ...
            'VariableNames',{'m2','even_w_min','E0','<x2>_even','<x4>_even','B_even', ...
            'odd_w_min','E1','<x2>_odd','<x4>_odd','B_odd'});
        writetable(outTable,[figdir,'mpx_d',num2str(d),'_delta',num2str(delta(ll)),'_lambda',num2str(lambda),'_cutoff',num2str(cutoff),'.txt'],'Delimiter','\t');
    end

    %% Finish plot
    xlabel('m2');
    ylabel('bc');

    xlim([m2(1),m2(end)]);

    cutoffLeg=repmat(cutoffAll,2,1);
    cutoffLeg=cutoffLeg(:);
    legend(string(cutoffLeg));

    title(['d: ',num2str(d),', delta: ',num2str(delta(ll)),', lambda: ',num2str(lambda)]);

    grid on
    box on

    print([figdir,'mpxFig_d',num2str(d),'_delta',num2str(delta(ll)),'_lambda',num2str(lambda),'.png'],'-dpng','-r0');

end

elapsedTime=toc;

disp(['Processing took ',num2str(elapsedTime/60/60,2),' hours.']);

%% Stop parallel pool

poolobj=gcp('nocreate');
delete(poolobj);