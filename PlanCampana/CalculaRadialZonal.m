Limpia
Region='CanaryIslands';

%Calculos en longitud
lon1=-15.3690;lat1=28.0830;
lon0=-14.5222;lat0=28.0830;

m=(lat1-lat0)/(lon1-lon0);

lonf=[-15.3690  -15.3249  -15.3008  -15.1888  -15.0537  -14.9185  -14.8045  -14.7104  -14.5800];
lonf=fliplr(lonf);
latf=lonf.*m-lon0*m+lat0;
for ii=1:length(latf)
       fprintf('%03d,%7.4f,%7.4f,1\n ',65+ii,lonf(ii),latf(ii))
end

%Configuracion geografica
lat_min=min(latf)-1; lat_max=max(latf)+1;  lon_min=360+min(lonf)-1;  lon_max=360+max(lonf)+1;
load(strcat(GlobalSU.LibPath,'/Settings/DS',Region)); %Cargo Data Settings
load(GlobalDS.filebat)

%Profundidad en las estacioens
for ii=1:length(latf)
    proe(ii)=-elevations(Locate(batylat,latf(ii)),Locate(batylon,lonf(ii)+360));
end

%Figuras
figure
plot(lonf,proe,'-o')

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
m_plot(360+lon0,lat0,'sr','markersize',16);
m_plot(360+lonf,latf,'ob-');

