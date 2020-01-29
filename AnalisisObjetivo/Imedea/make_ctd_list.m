function make_ctd_list(inpath,start,ending,outpath,...
   suf,temp,sal,sgth,pres,lat,lon,date);
%Este programa pasa los archivos de datos matlab a 
%formato *.prs, formato que será leido por los 
%programas fortran
%
%Eugenio Fraile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cambio los nan por 999 en todas las variables
ind=isnan(temp);temp(ind)=99.99;
ind=isnan(sal);sal(ind)=99.99;
ind=isnan(sgth);sgth(ind)=99.99;

outname=[outpath,'\info\ctd.lis'];
fid=fopen(outname,'w');

fprintf(fid,'######  LIST OF CTD_FILES USED TO CONSTRUCT ISOBARIC_LEVEL DATA FILES.  ######\n');
fprintf(fid,'######    INTERPOLATED_DATA FILES ARE (*.PRS) AND CONTAIN A HEADER      ######\n');
fprintf(fid,'\n');
for ii=1:size(lat,2)
   fprintf(fid,[suf,'d',sprintf('%3.3i',ii),'\n']);
end

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%genera los archivos *.prs en formato ascii para que
%puedan ser leidos por el programa fortran
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=1:size(lat,2)
   fl_prs=[outpath,'\prs\',suf,'d',sprintf('%3.3i',ii),'.prs'];
   fid=fopen(fl_prs,'w');
   
   %calculo lat en grados minutos y sg
   lat_g=floor(lat(ii));
   lat_m=floor((lat(ii)-lat_g)*60);
   lat_sg=(((lat(ii)-lat_g)*60)-lat_m)*60;
   
   %calculo lon en grados minutos y sg
   lon_g=fix(lon(ii));
   lon_m=-fix((lon(ii)-lon_g)*60);
   lon_sg=(-((lon(ii)-lon_g)*60)-lon_m)*60;
   
   fprintf(fid,[' St.Lat: ',num2str(lat_g),' 'sprintf('%2.0f',lat_m),...
         ' ',sprintf('%2.2f',lat_sg),' St.Lon:',num2str(lon_g),' ',...
			sprintf('%2.0f',lon_m),' ',sprintf('%2.2f',lon_sg),'\n']);
   
   fprintf(fid,'\n');
   fprintf(fid,'\n');
   fprintf(fid,'PRESS* TEMP* SALT* OXIG* SGTH* TRAN* FLU* NOBS\n');
   for jj=1:size(temp,1)
      fprintf(fid,[sprintf('%4.1f',pres(jj,ii)),' ',sprintf('%2.4f',temp(jj,ii)),' ',...
            sprintf('%2.4f',sal(jj,ii)),' ',sprintf('%2.2f',99.99),' ',...
            sprintf('%2.4f',sgth(jj,ii)),' ',sprintf('%2.2f',99.99),' ',...
            sprintf('%2.2f',99.99),' ',sprintf('%2.2f',99.99),'\n']);
   end
     
   fclose(fid);
end
