clear all
close all

tic

figdir='/scr/virga1/rsfdata/projects/qft/d5_deltaVaries_lambda1/';
%figdir='/scr/virga1/rsfdata/projects/qft/d/';

% Initialize
d=8;
nMax=10;

w=1:0.25:10;
lambda=1;
delta=1/d;
%delta=0.05:0.05:1;
m2=-(fliplr(0:2:20));
%m2=-(fliplr(12:1:16));
%m2=-14;

%cutoffAll=100:100:1000;
%cutoffAll=[100,500,1000,1100,1200,1300,1400,1500];
%cutoffAll=[5000,10000,15000];
cutoffAll=[20000,25000];
%cutoffAll=10000;

distSize=floor(d/2);

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
%uSorted=sort(C,2);
%uC=unique(uSorted,'rows');

sumC=sum(C,2);
uSumC=cat(2,C,sumC);

sC=sortrows(uSumC,size(uSumC,2));

close all
f1 = figure('Position',[200 500 600 600],'DefaultAxesFontSize',12,'renderer','painters');

t = tiledlayout(1,1,'TileSpacing','tight','Padding','tight');
s1=nexttile(1);
colmap=turbo(length(cutoffAll));

hold on

%% Loop over cutoff

bcE=nan(length(cutoffAll),length(m2));
bcO=nan(length(cutoffAll),length(m2));

cutoffLeg=[];

for kk=1:length(cutoffAll)
    cutoff=cutoffAll(kk);

    % Split into odd and even
    evensS=[];
    oddsS=[];

    for ii=1:cutoff
        if mod(sC(ii,end),2)==0
            evensS=cat(1,evensS,sC(ii,:));
        else
            oddsS=cat(1,oddsS,sC(ii,:));
        end
    end

    disp(['Cutoff: ',num2str(cutoff)]);

    disp('Even matrix');
    [matOutE1,matOutE2,matOutE3]=calcMatrix(evensS,d,w',delta);
    corrVecE=calcMatrix_corr(evensS,d);
    disp('Odd matrix');
    [matOutO1,matOutO2,matOutO3]=calcMatrix(oddsS,d,w',delta);
    corrVecO=calcMatrix_corr(oddsS,d);

    cutoffE=size(matOutE1,1);
    cutoffO=size(matOutO1,1);

    cutoffLeg=cat(1,cutoffLeg,cutoffE,cutoffO);

    minAllE=nan(1,length(m2));
    minAllEw=nan(1,length(m2));
    vecAllE=nan(cutoffE,length(m2));
    x2exE=nan(1,length(m2));
    x4exE=nan(1,length(m2));
    xCorrE=nan(length(m2),distSize);

    minAllO=nan(1,length(m2));
    minAllOw=nan(1,length(m2));
    vecAllO=nan(cutoffO,length(m2));
    x2exO=nan(1,length(m2));
    x4exO=nan(1,length(m2));
    xCorrO=nan(length(m2),distSize);

    for jj=1:length(m2)

        disp(['Loop over m: ',num2str(jj),' of ',num2str(length(m2))]);

        % Evens
        matOutE=matOutE1+matOutE2.*(m2(jj)+2./delta.^2)./2+matOutE3.*lambda./delta;

        % Eigen values
        minE=nan(1,length(w));
        vE=nan(cutoffE,length(w));

        for ii=1:length(w)
            %[V,D]=eig(matOutE(:,:,ii));
            [V,D]=eigs(matOutE(:,:,ii),2,'smallestabs');
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
        for bb=1:distSize
            xCorrE(jj,bb)=vecAllE(:,jj)'*corrVecE(:,:,bb)*vecAllE(:,jj)/(2*minAllEw(jj));
        end

        % Odds
        matOutO=matOutO1+matOutO2.*(m2(jj)+2./delta.^2)./2+matOutO3.*lambda./delta;

        % Eigen values
        minO=nan(1,length(w));
        vO=nan(cutoffO,length(w));

        for ii=1:length(w)
            %[V,D]=eig(matOutO(:,:,ii));
            [V,D]=eigs(matOutO(:,:,ii),2,'smallestabs');
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
        for bb=1:distSize
            xCorrO(jj,bb)=vecAllO(:,jj)'*corrVecO(:,:,bb)*vecAllO(:,jj)/(2*minAllOw(jj));
        end
    end

    bcE(kk,:)=1-x4exE./x2exE.^2./(3.*delta.*d);
    bcO(kk,:)=1-x4exO./x2exO.^2./(3.*delta.*d);

    plot(m2,bcE(kk,:),'-s','LineWidth',2,'Color',colmap(kk,:));
    plot(m2,bcO(kk,:),'--o','LineWidth',2,'Color',colmap(kk,:));
    drawnow

    outTable=table(m2',minAllEw',minAllE',x2exE',xCorrE,x4exE',bcE(kk,:)',minAllOw',minAllO',x2exO',xCorrO,x4exO',bcO(kk,:)', ...
        'VariableNames',{'m2','even_w_min','E0','<x2>_even','xCorr_even','<x4>_even','B_even', ...
        'odd_w_min','E1','<x2>_odd','xCorr_odd','<x4>_odd','B_odd'});
    writetable(outTable,[figdir,'mpx_d',num2str(d),'_delta',num2str(delta),'_lambda',num2str(lambda),'_cutoffTot',num2str(cutoff),'.txt'],'Delimiter','\t');
end

%% Finish plot
xlabel('m2');
ylabel('bc');

xlim([m2(1)-0.1,m2(end)+0.1]);

legend(string(cutoffLeg));

title(['d: ',num2str(d),', delta: ',num2str(delta),', lambda: ',num2str(lambda)]);

grid on
box on

print([figdir,'mpxFig_d',num2str(d),'_delta',num2str(delta),'_lambda',num2str(lambda),'_maxCutoof',num2str(cutoffAll(end)),'.png'],'-dpng','-r0');

elapsedTime=toc;

disp(['Processing took ',num2str(elapsedTime/60/60,2),' hours.']);

%% Stop parallel pool

poolobj=gcp('nocreate');
delete(poolobj);