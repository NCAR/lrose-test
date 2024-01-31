% Statistics of ecco output

clear all;
close all;

addpath(genpath('~/git/lrose-test/convstrat/dataProcessing/'));

startTime=datetime(2022,5,2,0,0,0);
endTime=datetime(2022,5,7,0,0,0);

saveData=1;

indir=['/scr/cirrus2/rsfdata/projects/nexrad-mrms/ecco_conus/'];

figdir=['/scr/cirrus2/rsfdata/projects/nexrad-mrms/figures/eccoStats/'];

flist=makeFileList(indir,startTime,endTime,'20YYMMDDxhhmmss',1);

timeBase=datetime(1970,1,1);

lon=ncread(flist{1},'x0');
lat=ncread(flist{1},'y0');
%alt=ncread(flist{1},'z2');

echoType2D.sl=int8(zeros(length(lat),length(lon),24));
echoType2D.sm=echoType2D.sl;
echoType2D.sh=echoType2D.sl;
echoType2D.m=echoType2D.sl;
echoType2D.ce=echoType2D.sl;
echoType2D.cs=echoType2D.sl;
echoType2D.cm=echoType2D.sl;
echoType2D.cd=echoType2D.sl;

%% Loop through radar files

for ii = 1:length(flist)
    file=flist{ii};

    disp(file);

    timeIn=ncread(file,'time');
    time=timeBase+seconds(timeIn);

    mst=time+minutes(lon.*4);
    mstHour=hour(mst);
    mstHour=repmat(mstHour',length(lat),1);

    echoType2Din=(ncread(file,'EchoTypeComp'))';

    % Loop through hours
    for jj=1:24
        addMat=echoType2Din==14 & mstHour==jj;
        echoType2D.sl(:,:,jj)=echoType2D.sl(:,:,jj)+int8(addMat);
    end

end