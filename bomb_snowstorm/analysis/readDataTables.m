function data=readDataTables(filename, del1)
% Reads Mikes txt tables and convert to matlab table

data.azimuth=[];
data.elevation=[];
data.range=[];
data.DBZ_F=[];
data.ZDR_F=[];
data.VEL_F=[];
data.PHIDP_F=[];
data.RHOHV_NNC_F=[];
data.WIDTH_F=[];
data.REGR_ORDER=[];
data.CMD_FLAG=[];
data.TRIP=[];

fid=fopen(filename,'r');
slurp=fscanf(fid,'%c');
fclose(fid);
first=1;

M=strread(slurp,'%s','delimiter','\n');

% Find first az= line
azInd=0;
searchInd=1;
while azInd==0
    thisLine=M{searchInd};
    if strcmp(thisLine(1:2),'az')
        azInd=searchInd;
    end
    searchInd=searchInd+1;
end

for ii=azInd:length(M)
    thisLine=M{ii};
    if strcmp(thisLine(1:2),'az')
        thisStr=strsplit(thisLine,del1);
        data.azimuth=cat(1,data.azimuth,str2double(thisStr{2}));
        data.elevation=cat(1,data.elevation,str2double(thisStr{4}));
        if first==0;
            data.range=cat(1,data.range,rayMat(:,1)');
            data.DBZ_F=cat(1,data.DBZ_F,rayMat(:,2)');
            data.ZDR_F=cat(1,data.ZDR_F,rayMat(:,3)');
            data.VEL_F=cat(1,data.VEL_F,rayMat(:,4)');
            data.PHIDP_F=cat(1,data.PHIDP_F,rayMat(:,5)');
            data.RHOHV_NNC_F=cat(1,data.RHOHV_NNC_F,rayMat(:,6)');
            data.WIDTH_F=cat(1,data.WIDTH_F,rayMat(:,7)');
            data.REGR_ORDER=cat(1,data.REGR_ORDER,rayMat(:,8)');
            data.CMD_FLAG=cat(1,data.CMD_FLAG,rayMat(:,9)');
            if size(rayMat,2)==10
                data.TRIP=cat(1,data.TRIP,rayMat(:,10)');
            end
        end
        rayMat=[];
    else
        temp=strread(M{ii},'%f','delimiter',del1);
        rayMat=cat(1,rayMat,temp');
        first=0;
    end
end
data.range=cat(1,data.range,rayMat(:,1)');
data.DBZ_F=cat(1,data.DBZ_F,rayMat(:,2)');
data.ZDR_F=cat(1,data.ZDR_F,rayMat(:,3)');
data.VEL_F=cat(1,data.VEL_F,rayMat(:,4)');
data.PHIDP_F=cat(1,data.PHIDP_F,rayMat(:,5)');
data.RHOHV_NNC_F=cat(1,data.RHOHV_NNC_F,rayMat(:,6)');
data.WIDTH_F=cat(1,data.WIDTH_F,rayMat(:,7)');
data.REGR_ORDER=cat(1,data.REGR_ORDER,rayMat(:,8)');
data.CMD_FLAG=cat(1,data.CMD_FLAG,rayMat(:,9)');
if size(rayMat,2)==10
    data.TRIP=cat(1,data.TRIP,rayMat(:,10)');
else
    data=rmfield(data,'TRIP');
end

data.range=data.range(1,:);
end