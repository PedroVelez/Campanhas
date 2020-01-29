Limpia
Region='CanaryIslands';

% %Calculos en latitud
lon0=-14.0000;lat0=26.5000;
lon1=-15.5706; lat1=27.7053;
m=(lat1-lat0)/(lon1-lon0);
latf=[27.7053 27.6574 27.5836 27.4828 27.3546   27.2275   27.1059 26.9476 26.8724 26.8009 26.7293 26.6691 26.5861];
latf=fliplr(latf);

% %Calculos en latitud
% lon1=-14.0000;lat1=28.1463;
% lon0=-14.0000;lat0=26.5000;
% m=(lat1-lat0)/(lon1-lon0);
% latf=[28.1700 28.1185 28.0700 28.0033   27.8083   27.6133   27.4183   27.2833 27.14 27.0283 26.9476 26.8724 26.8009 26.7293 26.6691 26.5861 26.5000];

lonf=lon0+(latf-lat0)/m;
for ii=1:length(latf)
    fprintf('%03d,%7.4f,%7.4f,1\n',34+ii,lonf(ii),latf(ii))
end

%Configuracion geografica
lat_min=min(latf)-0.5; lat_max=max(latf)+0.5;  lon_min=360+min(lonf)-0.5;  lon_max=360+max(lonf)+0.5;
load(strcat(GlobalSU.LibPath,'/Settings/DS',Region)); %Cargo Data Settings
load(GlobalDS.filebat)

%Profundidad en las estacioens
for ii=1:length(latf)
    proe(ii)=-elevations(Locate(batylat,latf(ii)),Locate(batylon,lonf(ii)+360));
end

%Figuras
figure
plot(latf,proe,'-o')

figure
m_proj('mercator','long',[lon_min lon_max],'lat',[lat_min lat_max]);hold on
m_contour(batylon,batylat,elevations,[-3000 -3000],'color',[0.45 0.45 0.45]);
m_contour(batylon,batylat,elevations,[-2000 -2000],'color',[0.55 0.55 0.45]);
m_contour(batylon,batylat,elevations,[-1000 -1000],'color',[0.65 0.65 0.65]);
m_contour(batylon,batylat,elevations,[ -500  -500],'color',[0.75 0.75 0.75]);
m_contour(batylon,batylat,elevations,[ -250  -250],'color',[0.85 0.85 0.85]);
m_usercoast(GlobalDS.filecoast,'patch',[.7 .6 .4,],'edgecolor',[.7 .6 .4,]);
m_grid

m_plot(360+lon1,lat1,'sr','markersize',6);hold on
m_plot(360+lon0,lat0,'sr','markersize',6);
m_plot(360+lonf,latf,'ob-');

