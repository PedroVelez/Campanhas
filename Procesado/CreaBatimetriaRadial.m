Limpia

%FileData='Ra1804LaGraciosa';Opciones.tipo='m';
%FileData='Ra1804CaboGhir';Opciones.tipo='z';
%FileData='Ra1804Norte';Opciones.tipo='z';
%FileData='Ra1804LP';Opciones.tipo='m';
%FileData='Ra1804Agadir';Opciones.tipo='z';
FileData='Ra1804CJuby';Opciones.tipo='z';

FileOut=strcat('Batimetria',FileData,'.mat');
FileBat='CanaryIslandsHRBat';
Opciones.disnt=.1; %Distancia [km] puntos a crear entre cada estacion en base a la batimetria

%%
if strcmp(Opciones.tipo,'Zonal')
    Opciones.tipo='Z';
    fprintf('         > Zonal section\n');
elseif strcmp(Opciones.tipo,'zonal')
    Opciones.tipo='Z';
    fprintf('         > Zonal section\n');
elseif strcmp(Opciones.tipo,'z')
    Opciones.tipo='Z';
    fprintf('         > Zonal section\n');
elseif strcmp(Opciones.tipo,'Meridional')
    Opciones.tipo='M';
    fprintf('         > Meridional section\n');
elseif strcmp(Opciones.tipo,'meridional')
    Opciones.tipo='M';
    fprintf('         > Meridional section\n');
elseif strcmp(Opciones.tipo,'m')
    Opciones.tipo='M';
    fprintf('         > Meridional section\n');
else
    Opciones.tipo='Z';
    fprintf('         > Zonal section [Default]\n');
end
if ~isfield(Opciones,'tipo')
    Opciones.tipo='Z';
    fprintf('         > Zonal section [Default]\n');
end


load(FileBat)
DATA=load(FileData);
IndSec=1:length(DATA.longs);

lon=DATA.longs(IndSec)+360;
lat=DATA.latis(IndSec);

for i1=1:length(lon)
    ilat1=Locate(batylat,lat(i1));
    ilon1=Locate(batylon,lon(i1));
    pro(i1)=elevations(ilat1,ilon1);
end
latit=[];
lonit=[];
%Construyo una vector con mayor definicion que el vector de posiciones de las estaciones
for i=2:length(lon)
    dist=sw_dist([lat(i-1) lat(i)],[lon(i-1) lon(i)],'km');
    lonit=[lonit interp1([1 2],[lon(i-1) lon(i)],linspace(1,2,dist/Opciones.disnt))];
    latit=[latit interp1([1 2],[lat(i-1) lat(i)],linspace(1,2,dist/Opciones.disnt))];
end
for ii=1:length(latit)
    ilat1=Locate(batylat,latit(ii));
    ilon1=Locate(batylon,lonit(ii));
    x=[batylat(ilat1) batylat(ilat1+1); batylat(ilat1) batylat(ilat1+1)];
    y=[batylon(ilon1) batylon(ilon1);   batylon(ilon1+1) batylon(ilon1+1)];
    z=[elevations(ilat1,ilon1) elevations(ilat1+1,ilon1);
        elevations(ilat1,ilon1+1) elevations(ilat1+1,ilon1+1)];
    zi(ii)=interp2(x,y,z,latit(ii),lonit(ii));
end



batylat=latit';
batylon=lonit';
elevations=zi;

figure
load('../DatosCampanha')
m_proj('Mercator','long',[lon_min lon_max],'lat',[lat_min lat_max]);
m_plot(lon,lat,'ob');hold on
m_plot(batylon,batylat,'r.');hold on
if ~isempty(filecosta)
    m_usercoast(filecosta,'patch',[.7 .6 .4,],'edgecolor',[.7 .6 .4,]);
end
m_grid('linestyle','none','fontsize',08,'fontname','tahona')

figure
if  Opciones.tipo=='M';
    area(batylat,elevations)
    hold on;plot(lat,pro,'s','markersize',12)
else
    area(batylon-360,elevations)
    hold on;plot(lon-360,pro,'s','markersize',12)
end


save(FileOut,'elevations','batylon','batylat')