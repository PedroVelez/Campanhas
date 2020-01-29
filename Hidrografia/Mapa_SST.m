close all;clear all;clc

DC=load('../DatosCampanha');

SVOpciones.DataSet=DC.campanhacode;

SVOpciones.SateliteDataSet='noaa.oisst.v2.highres';
SVOpciones.datei=datenum(2018,04,06);
SVOpciones.Caxis=0;[14.5 22];
SVOpciones.CaxisAnomaly=[-3 1];
SVOpciones.figuras=[4 7];
SVOpciones.colormap='jet';
SVOpciones.SatellitePath=strcat('/Users/pvb/Dropbox/Oceanografia/Data/Satelite/',SVOpciones.SateliteDataSet,'/');
SVOpciones.SatelliteFile=strcat(SVOpciones.SatellitePath,'NOAAOisstv2HighresSstDayMean_01011982_06042018_CCLME.mat');

SVOpciones.filecosta=DC.filecosta;
SVOpciones.lon_min=DC.lon_minEXT;
SVOpciones.lon_max=DC.lon_maxEXT;
SVOpciones.lat_min=DC.lat_minEXT;
SVOpciones.lat_max=DC.lat_maxEXT;


%% Read data
SVOpciones.ztitulo=SVOpciones.SateliteDataSet;
SVOpciones.titulo=DC.campanha;

[Ye,Mo,Da]=datevec(SVOpciones.datei);
Data=matfile(SVOpciones.SatelliteFile);

ilon_min=Locate(Data.lon,SVOpciones.lon_min);
ilon_max=Locate(Data.lon,SVOpciones.lon_max);
ilat_min=Locate(Data.lat,SVOpciones.lat_min);
ilat_max=Locate(Data.lat,SVOpciones.lat_max);
ijday=Locate(Data.timetd,SVOpciones.datei);

lon=Data.lon(ilon_min:ilon_max,1);
lat=Data.lat(ilat_min:ilat_max,1);
jday=Data.timetd(ijday,1);
jdayT=Data.timetd(:,1);
[YT,MT,DT]=datevec(jdayT);

zvarT=Data.ssttd(ilon_min:ilon_max,ilat_min:ilat_max,:);
zvar=Data.ssttd(ilon_min:ilon_max,ilat_min:ilat_max,ijday);

zvarm=nanmean(zvarT(:,:,MT==Mo),3);


%% figure datei
figure
fprintf('     >%s %s, Minx:%6.3f, Max:%6.3f\n',SVOpciones.titulo,SVOpciones.ztitulo,min(zvar(:)),max(zvar(:)))
m_proj('mercator','lon',[SVOpciones.lon_min SVOpciones.lon_max],'lat',[SVOpciones.lat_min SVOpciones.lat_max])
[c,h1]=m_contourf(lon,lat,zvar',80,'edgecolor','none');hold on
XL=nanmin(lon(:))+0.05*Rango(lon);
YL=nanmax(lat(:))-0.025*Rango(lat);
h3=m_text(XL,YL,sprintf('%4d %02d %02d',Ye,Mo,Da),'backgroundcolor','w','Fontsize',16);
%estaciones
DataCrusie=load(SVOpciones.DataSet);
for i1=1:size(DataCrusie.longs,2)
    m_plot(DataCrusie.longs(i1)+360,DataCrusie.latis(i1),'marker','o','markersize',4,'color','k','MarkerFaceColor','w');
    
end
title(sprintf('%s %s',SVOpciones.titulo,SVOpciones.ztitulo),'FontSize',12,'Fontweight','bold','interpreter','none');
m_usercoast(DC.filecosta,'patch',[.7 .6 .4,],'edgecolor',[.7 .6 .4,]);hold on
m_grid('linestyle','none')
colorbar
if SVOpciones.Caxis==0
    colormap(SVOpciones.colormap);
else
    caxis([SVOpciones.Caxis])
    colormap(SVOpciones.colormap);
end

orient landscape;CreaFigura(gcf,deblankT(strcat(SVOpciones.titulo,SVOpciones.ztitulo,datestr(SVOpciones.datei,'YYYYmmDD'))),SVOpciones.figuras)


%% promedio mensual
figure
fprintf('     >%s %s Monthly Minx:%6.3f, Max:%6.3f\n',SVOpciones.titulo,SVOpciones.ztitulo,min(zvarm(:)),max(zvarm(:)))
m_proj('mercator','lon',[SVOpciones.lon_min SVOpciones.lon_max],'lat',[SVOpciones.lat_min SVOpciones.lat_max])
[c,h1]=m_contourf(lon,lat,zvarm',80,'edgecolor','none');hold on
title(sprintf('Monthly (%d) mean %s %s',Mo,SVOpciones.titulo,SVOpciones.ztitulo),'FontSize',12,'Fontweight','bold','interpreter','none');
m_usercoast(DC.filecosta,'patch',[.7 .6 .4,],'edgecolor',[.7 .6 .4,]);hold on
m_grid('linestyle','none')
colorbar

if SVOpciones.Caxis==0
    colormap(SVOpciones.colormap);
else
    caxis([SVOpciones.Caxis])
    colormap(SVOpciones.colormap);
end
orient landscape;CreaFigura(gcf,deblankT(strcat(SVOpciones.titulo,SVOpciones.ztitulo,'MonthlyMean')),SVOpciones.figuras)

%% anomaly
figure
fprintf('     >%s %s Anomaly (SST-Monthly mean), Minx:%6.3f, Max:%6.3f\n',SVOpciones.titulo,SVOpciones.ztitulo,min(zvar(:)-zvarm(:)),max(zvar(:)-zvarm(:)))
m_proj('mercator','lon',[SVOpciones.lon_min SVOpciones.lon_max],'lat',[SVOpciones.lat_min SVOpciones.lat_max])
[c,h1]=m_contourf(lon,lat,zvar'-zvarm',80,'edgecolor','none');hold on
titulo=title(sprintf('Anomaly (SST-Monthly(%02d)mean)  %s %s',Mo,SVOpciones.titulo,SVOpciones.ztitulo),'FontSize',12,'Fontweight','bold','interpreter','none');
m_usercoast(DC.filecosta,'patch',[.7 .6 .4,],'edgecolor',[.7 .6 .4,]);hold on
m_grid('linestyle','none')
colorbar
caxis([SVOpciones.CaxisAnomaly])
colormap(SVOpciones.colormap);
orient landscape;CreaFigura(gcf,deblankT(strcat(SVOpciones.titulo,SVOpciones.ztitulo,datestr(SVOpciones.datei,'YYYYmmDD'),'Anomaly')),SVOpciones.figuras)

