function [lat,lon,nstat,gtime,depth,data,names]=Sbe911_mat(file);
%
% SBE911_mat Reads the SeaBird ASCII .cnv
%
%     Usage: [lat,lon,nstat,ncast,gtime,depth,data,names,sensors]=sbe911_mat(file);
%
%     Input:  file = name of file  (e.g. 'cast002.cnv')
%
%     Output: lon = longitude in decimal degrees, West negative
%             lat = latitude in decimal degrees, North positive
%           gtime = Gregorian time vector in UTC
%           depth = Bottom depth
%            data = matrix containing all the columns of data in the .CNV file
%           names = string matrix containing the names and units of the columns
%
% Open the .cnv file as read-only text
%
fid=fopen(file,'rt');
%
% Read the header.
% Start reading header lines of .CNV file,
% Stop at line that starts with '*END*'
%
% Pull out lat & lon along the way and look
% at the '# name' fields to see how many variables we have.
%
str='*START*'; %#ok<NASGU>
var=0;
sens=0;
str=fgetl(fid);

while (~strncmp(str,'*END*',5));
    str=fgetl(fid);
    %-----------------------------------
    if ~isempty(strfind(str,'Station:'))
        is=strfind(str,':');
        nstat=str2double(str(is+1:end));
        %     elseif ~isempty(strfind(str,'Latitude:'))
        %         is=min(strfind(str,':'));
        %         dm=sscanf(str(is+1:end),'%d:%f %1s');
        %         lat=dm(1)+(dm(2))/60;
        %     elseif ~isempty(strfind(str,'Longitude:'))
        %         is=min(strfind(str,':'));
        %         dm=sscanf(str(is+1:end),'%d:%f %s');
        %         lon=(abs(dm(1))+dm(2)/60)*sign(dm(1));
        %     elseif (strfind(str,'* System UpLoad')==1)
        %         is=strfind(str,'=');
        %         datstr=[str(is+6:is+7) '-' str(is+2:is+4) '-' str(is+9:is+12)];
        %         n=datenum(datstr);
        %         gtime=datevec(n);
        %         gtime(4)=str2double(str(is+14:is+15));
        %         gtime(5)=str2double(str(is+17:is+18));
        %         gtime(6)=str2double(str(is+20:is+21));
    elseif (strfind(str,'* NMEA Latitude')==1)
        is=strfind(str,'=');
        dm=sscanf(str(is+1:end),'%d %f %1s');
        latNMEA=dm(1)+(dm(2))/60;
        if strcmp(char(dm(3)),'S') || strcmp(char(dm(3)),'s')
            latNMEA=latNMEA;
        end
    elseif (strfind(str,'* NMEA Longitude')==1)
        is=strfind(str,'=');
        dm=sscanf(str(is+1:end),'%d %f %s');
        lonNMEA=(abs(dm(1))+dm(2)/60);
        if strcmp(char(dm(3)),'W')|| strcmp(char(dm(3)),'w')
            lonNMEA=-lonNMEA;
        end

    elseif (strfind(str,'* NMEA UTC (Time) =')==1)
        is=strfind(str,'=');
        datstr=[str(is+6:is+7) '-' str(is+2:is+4) '-' str(is+9:is+12)];
        n=datenum(datstr);
        gtimeNMEA=datevec(n);
        iss=strfind(str,':');
        gtimeNMEA(4)=str2num(str(iss(1)-2:iss(1)-1));
        gtimeNMEA(5)=str2num(str(iss(2)-2:iss(2)-1));
        gtimeNMEA(6)=str2num(str(iss(2)+1:iss(2)+2));
    elseif ~isempty(strfind(str,'** Depth'))
        is=strfind(str,':');
        dm=sscanf(str(is+1:end),'%d,%d');
        if isempty(dm)
            depth=NaN;
            disp('     There is not Depth values in the cnv file')
        else
            depth=dm(1);
        end
    elseif ~isempty(strfind(str,'Bottom Depth [m]:'))
        is=strfind(str,':');
        dm=sscanf(str(is+1:end),'%d,%d');
        depth=dm(1);
    elseif (strfind(str,'# name')==1)
        var=var+1;
        names{var}=str;
    elseif (strfind(str,'# bad_flag')==1)
        isub=13:length(str);
        bad_flag=sscanf(str(isub),'%g',1);
    end
end
if exist('lonNMEA','var')==1
    lon=lonNMEA;
    disp('     Using NMEA lon')
end
if exist('latNMEA','var')==1
    lat=latNMEA;
    disp('     Using NMEA lat')
end
if exist('gtimeNMEA','var')==1
    gtime=gtimeNMEA;
    disp('     Using NMEA gtime')
end
%==============================================
%  Done reading header.  Now read the data!
%==============================================
e
nvars=var;  %number of variables
% Read the data into one big matrix
data=fscanf(fid,'%f',[nvars inf]);
fclose(fid);

% Flag bad values with nan
ind=find(data==bad_flag);
data(ind)=data(ind)*nan;

% Flip data around so that each variable is a column
data=data.';
% Convert cell arrays of names to character matrices
names=char(names);

return
