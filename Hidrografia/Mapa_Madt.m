close all;clear all;clc

DC=load('../DatosCampanha');

SVOpciones.DataSet=DC.campanhacode;

SVOpciones.SateliteDataSet='nrt_global_allsat_phy_l4';
SVOpciones.datei=datenum(2018,04,02);
SVOpciones.facm=0.02;
SVOpciones.Caxis=[20 48];
SVOpciones.figuras=[4 7];
SVOpciones.colormap='jet';
SVOpciones.SatellitePath=strcat('/Users/pvb/Dropbox/Oceanografia/Data/Satelite/Aviso/',SVOpciones.SateliteDataSet,'/Mat/');

SVOpciones.filecosta=DC.filecosta;
SVOpciones.lon_min=DC.lon_minEXT;
SVOpciones.lon_max=DC.lon_maxEXT;
SVOpciones.lat_min=DC.lat_minEXT;
SVOpciones.lat_max=DC.lat_maxEXT;

%% Read data
SVOpciones.ztitulo=SVOpciones.SateliteDataSet;
SVOpciones.titulo=DC.campanha;

[Ye,Mo,Da]=datevec(SVOpciones.datei);

Data=load(strcat(SVOpciones.SatellitePath,'/',SVOpciones.SateliteDataSet,sprintf('_%4d%02d%02d',Ye,Mo,Da)))
X=double(Data.X);
Y=double(Data.Y);

%Recorte la zona
ilon_min=Locate(X(1,:),SVOpciones.lon_min);
ilon_max=Locate(X(1,:),SVOpciones.lon_max);
ilat_min=Locate(Y(:,1),SVOpciones.lat_min);
ilat_max=Locate(Y(:,1),SVOpciones.lat_max);

lon(:,:)=X(ilat_min:ilat_max,ilon_min:ilon_max);
lat(:,:)=Y(ilat_min:ilat_max,ilon_min:ilon_max);
jday=Data.julianday;

zvar(:,:)=Data.madt(ilat_min:ilat_max,ilon_min:ilon_max)*100;
ur(:,:)=Data.u(ilat_min:ilat_max,ilon_min:ilon_max)*100;
vr(:,:)=Data.v(ilat_min:ilat_max,ilon_min:ilon_max)*100;

%% figure
figure
m_proj('mercator','lon',[DC.lon_minEXT DC.lon_maxEXT],'lat',[DC.lat_minEXT DC.lat_maxEXT])
[c,h1]=m_contourf(lon,lat,zvar(:,:),40,'edgecolor','none');hold on
h2=m_quiver(lon,lat,ur(:,:)*SVOpciones.facm,vr(:,:)*SVOpciones.facm,0,'filled','k');hold on
XL=nanmin(lon(:))+0.05*Rango(lon);
YL=nanmax(lat(:))-0.025*Rango(lat);
h3=m_text(XL,YL,sprintf('%4d %02d %02d',Ye,Mo,Da),'backgroundcolor','w','Fontsize',16);

%estaciones
DataCrusie=load(SVOpciones.DataSet);
for i1=1:size(DataCrusie.longs,2)
    m_plot(DataCrusie.longs(i1)+360,DataCrusie.latis(i1),'marker','o','markersize',8,'MarkerFaceColor','w');
    
end
title(sprintf('%s %s',SVOpciones.titulo,SVOpciones.ztitulo),'FontSize',12,'Fontweight','bold','interpreter','none');
m_usercoast(DC.filecosta,'patch',[.7 .6 .4,],'edgecolor',[.7 .6 .4,]);hold on
m_grid('linestyle','none')
colorbar
caxis([SVOpciones.Caxis])
colormap(SVOpciones.colormap);

fprintf('     >%s %s, Max:%6.3f, Min:%6.3f\n',SVOpciones.titulo,SVOpciones.ztitulo,max(zvar(:)),min(zvar(:)))

orient landscape;CreaFigura(gcf,deblankT(strcat(SVOpciones.titulo,SVOpciones.ztitulo,datestr(SVOpciones.datei,'YYYYmmDD'))),SVOpciones.figuras)
