%Script para el post-procesamiento de los datos obtenidos
%del análisis objetivo
%
%Eugenio Fraile
%efraile@becarios.ulpgc.es
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%ruta a la carpeta analisis_objetivo
outpath='D:\misdocumentos\Eugenio\Gyroscope\gyros0302\analisis_objetivo2\';

%cargo los datos iniciales del archivo grid.mat
eval(['cd ',outpath,'info\']);
load grid

%Esta función extrae de toda la malla de puntos, las
%estaciones de la periferia de la caja.
extrae_sta_peri(xlon1,xlat1,NL,NC,NF,...
   xarm,yarm,zo,zint,ztot,suf,outpath);

