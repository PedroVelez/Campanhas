close all;clear all;clc
%---------------------------------------------------
%Guardo los datos para visualizacion en ODV
%---------------------------------------------------

load('DatosCampanha')

%--------------------------------
datadir=strcat(dirdata,'\datos\Mat\');
disp(sprintf('>>>> Pasando de formato Mat a ODV'))
[file,pathfile]=uigetfile ( '*.mat', 'Escoge la seccion a representar' );load(strcat(pathfile,file))

fileout=strcat(file(1:end-4),'_odv.dat');

if exist('oxygs')==0
    oxygs=temps.*NaN;
end
if exist('fluos')==0
    oxygs=fluos.*NaN;
end
[year,month,day,hour,min,sec]=datevec(dates);

temps(isnan(temps))=-999;
salts(isnan(salts))=-999;
sgths(isnan(sgths))=-999;
fluos(isnan(fluos))=-999;
oxygs(isnan(oxygs))=-999;

fod=fopen(fileout,'w');
%Construyo el cabecero para el fichero de datos woce.
%fprintf(fod,'Station ;Year ;Month ;Day ;Latitude ;Longitude ;Pressure [db] ;Temperature [ITS-90, deg C] ;Salinity [PSU] ;Fluorescence ;Oxygen [ml/l] \n');
fprintf(fod,'Station;Year;Month;Day;Latitude;Longitude;Pressure[db];Temperature[ITS-90, deg C];Salinity;Pdensity[Kg/m3];Fluorescence,Seapoint;Oxygen[ml/l]\n');
for ii=1:size(nstats,2)
    y=ones(size(temps,1),1)*year(ii);
    m=ones(size(temps,1),1)*month(ii);
    d=ones(size(temps,1),1)*day(ii);
    st=ones(size(temps,1),1)*nstats(ii);
    lati=ones(size(temps,1),1)*latis(ii);
    long=ones(size(temps,1),1)*longs(ii);
    %Guarda los datos formateados en un fichero ascii.
    guardar=[st y m d lati long press(:) temps(:,ii) salts(:,ii) sgths(:,ii) oxygs(:,ii) fluos(:,ii)];
    fprintf(fod,'%03d;%04d;%02d;%02d;%4.3f;%4.3f;%7.3f;%7.3f;%7.3f;%7.3f;%7.3f;%7.3f\n',guardar');
end
fclose(fod);
disp(sprintf('   > Grabando en %s',fileout))


