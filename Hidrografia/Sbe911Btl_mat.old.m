function [lat,lon,nstat,gtime,depth,data,names]=Sbe911Btl_mat(file)

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
S=regexp(str,'Bottle        Date');
%while (~strncmp(str,'  Position',10));
while isempty(S);
    str=fgetl(fid);
    S=regexp(str,'Bottle        Date');
    
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
        lonNMEA=(abs(dm(1))+dm(2)/60)*sign(dm(1));
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
    elseif ~isempty(strfind(str,'Depth:'))
        is=strfind(str,':');
        dm=sscanf(str(is+1:end),'%d,%d');
        depth=dm(1);
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

names=str;
str2C=fgetl(fid);
% Read Data
ii=0;
eofstat=feof(fid);
while ~eofstat
    ii=ii+1;
    str1=fgetl(fid);
    str2=fgetl(fid);
    t=sscanf(str1(26:end-5),'%f');
    
    S1=regexp(str1,'\s+','split');
    S=regexp(str,'\s+','split');
    for is=2:length(S)
        if strcmp(S{is},'Bottle')
            iBot=is-1;
        elseif strcmp(S{is},'Sal00')
            iSal=is-1;
        elseif strcmp(S{is},'PrDM')
            iP=is-1;
        elseif strcmp(S{is},'T090C')
            iTem=is-1;
        elseif strcmp(S{is},'C0S/m')
            iCon=is-1;
        elseif strcmp(S{is},'C1S/mTurbWETntu0')
            iTur=is-1+1; %!!Cuidado que debo suar uno porque C1S y TurbWETntu0 aparecen pegadps
        elseif strcmp(S{is},'SeaTurbMtr')
            iTur2=is-1+1;
        elseif strcmp(S{is},'FlECO-AFL')
            iFlu=is-1+1;
        elseif strcmp(S{is},'Sbeox0ML/L')
            iOxy=is-1+1;
        elseif strcmp(S{is},'Sbox0Mm/Kg')
            iOxy2=is-1+1;
            
            
        end
    end
    data.numbtl(ii)=str2double(S1{2});
    data.salbtl(ii)=str2double(cell2mat(S1(iSal+3)));
    data.prebtl(ii)=str2double(cell2mat(S1(iP+3)));
    data.tembtl(ii)=str2double(cell2mat(S1(iTem+3)));
    data.conbtl(ii)=str2double(cell2mat(S1(iCon+3)));
    if exist('iOxy')
        data.oxybtl(ii)=str2double(cell2mat(S1(iOxy+3)));
    end
    if exist('iOxy2')
        data.oxy2btl(ii)=str2double(cell2mat(S1(iOxy2+3)));
    end
    if exist('iTur')
        data.turbtl(ii)=str2double(cell2mat(S1(iTur+3)));
    end
    if exist('iTur2')
        data.tur2btl(ii)=str2double(cell2mat(S1(iTur2+3)));
    end
    if exist('iFlu')
        data.flubtl(ii)=str2double(cell2mat(S1(iFlu+3)));
    end
    eofstat=feof(fid);
end %while
fclose(fid);

return