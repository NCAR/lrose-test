function data=readData_Torres(filename)
% Reads Johns Torres tables and convert to matlab table

fid=fopen(filename,'r');
slurp=fscanf(fid,'%c');
fclose(fid);
first=1;

M=strread(slurp,'%s','delimiter','\n');

% Find first 64 line
sfInd=0;
searchInd=1;
while sfInd==0
    thisLine=M{searchInd};
    if strcmp(thisLine(1:3),'64 ')
        sfInd=searchInd;
    end
    searchInd=searchInd+1;
end

data.width=[];
data.tableOut=[];

for ii=sfInd:length(M)
    thisLine=M{ii};
    if strcmp(thisLine(1:3),'64 ')
        thisStr=strsplit(thisLine,' ');
        data.width=cat(1,data.width,str2double(thisStr{16}));        
        if first==0
            data.tableOut=cat(3,data.tableOut,rayMat);
        end
        rayMat=[];
    else
        temp=strread(M{ii},'%f','delimiter',' ');
        rayMat=cat(1,rayMat,temp');
        first=0;
    end
end
data.tableOut=cat(3,data.tableOut,rayMat);
end