%		DatosCampanha
%		Fichero donde se almacena toda la informacion de la campanha en curso
%		Usado por todos los programas del paquete HidroIEO
%
% 1.0 26 Junio 2004 Instituto Espa?ol de Oceanografia (c)
close all;clear all;clc

campanha='RaProCan1804';			% Nombre de la campanha
campanhacode='Ra1804';            % Codigo de la campanha

dirdata='/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/2018_Raprocan1804/CTD';
filebat=strcat('CanaryIslandsBat.mat');  	% Fichero con la batimetria de la zona
filecosta=strcat('CanaryIslandsCoast');	% Fichero con la costa de la zona

%Area geografica
lat_min=27;   lat_max=32;
lon_min=360-19;  lon_max=360-8;
figure
m_proj('mercator','lon',[lon_min lon_max],'lat',[lat_min lat_max])
title(sprintf('%s %s',campanha,campanhacode),'FontSize',12,'Fontweight','bold','interpreter','none');
m_usercoast(filecosta,'patch',[.7 .6 .4,],'edgecolor',[.7 .6 .4,]);hold on
m_grid('linestyle','none')
CreaFigura(gcf,strcat('DatosCampanha',campanha),4)

%Area geografica Extendida (SST,MADT)
lat_minEXT=22;   lat_maxEXT=34;
lon_minEXT=360-24;  lon_maxEXT=360-8;
figure
m_proj('mercator','lon',[lon_minEXT lon_maxEXT],'lat',[lat_minEXT lat_maxEXT])
title(sprintf('Extended %s %s',campanha,campanhacode),'FontSize',12,'Fontweight','bold','interpreter','none');
m_usercoast(filecosta,'patch',[.7 .6 .4,],'edgecolor',[.7 .6 .4,]);hold on
m_grid('linestyle','none')

CreaFigura(gcf,strcat('DatosCampanha',campanha,'EXT'),4)


save('DatosCampanha')