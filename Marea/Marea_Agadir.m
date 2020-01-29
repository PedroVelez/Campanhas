clear all; close all;

stations=[1110 914:-1:901];

CruiseDir='/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/2018_Raprocan1804';
output_file=strcat(CruiseDir,'/Analisis/Marea/Raprocan1804_marea_Agadir');
datadir=    strcat(CruiseDir,'/LADCP/Visbeck/profiles/');

for ii=1:length(stations)
    if stations(ii)<99
        sta=sprintf('Ra1804_%03d',stations(ii));
    else
        sta=sprintf('Ra1804_%03d',stations(ii));
    end
    ext='.mat';
    filename=[datadir,sta,ext];
    load(filename);
    
    lat(ii)=dr.lat;
    lon(ii)=dr.lon;
    time(ii)=datenum(dr.date);
end

% [TS_u,ConList_1]=tide_pred('Model_tpxo6.2',time,lat,lon+360,'u');
% [TS_v,ConList_1]=tide_pred('Model_tpxo6.2',time,lat,lon+360,'v');

TS_u=[];
TS_v=[];
for ii=1:length(stations)
    ii
    [TS_u_l,ConList_1_l]=tide_pred('/Users/pvb/Dropbox/Oceanografia/Data/Mareas/TPXO7.2/Model_tpxo7.2',time(ii),lat(ii),lon(ii)+360,'u');
    [TS_v_l,ConList_1_l]=tide_pred('/Users/pvb/Dropbox/Oceanografia/Data/Mareas/TPXO7.2/Model_tpxo7.2',time(ii),lat(ii),lon(ii)+360,'v');
    TS_u=merge(TS_u,TS_u_l);
    TS_v=merge(TS_v,TS_v_l);
end

distance=[sw_dist(lat,lon,'km')];
dist=[0 cumsum(distance)];

figure
% plot(lon,TS_v,'*-','linewidth',2) %s?lo componente zonal, el transecto es meridional
plot(dist,TS_v,'*-','linewidth',2) %s?lo componente zonal, el transecto es meridional
% xlim([min(lon) max(lon)])
ylim([-10 +10])
set(gca,'TickDir','in')
set(gca,'fontsize',12)
% xlabel('Longitude ({\circ}W)','fontsize',14)
xlabel('Distance (km)','fontsize',14)
ylabel('Velocity (cm/s)','fontsize',14)
grid on

%    ad1=add_scale_isis2('t',dist,stations);
%    set(ad1,'visible','off');
%     set(ad1,'xdir','rev')
%     set(gca,'tickdir','out');
a=get(gcf,'child');
c=get(a(1),'child');
tx=findobj(c,'type','text');
%         for ii=[1:length(stations)]
%          set(tx(ii),'fontsize',11);
%         end

htitle=title('Tidal component predicted from TPXO 7.2 global model','fontsize',16);
posi=get(htitle,'pos');
set(htitle,'pos',[posi(1) posi(2)+0.5 posi(3)]);
orient landscape

CreaFigura(1,output_file,[7])
save(output_file,'TS_u','TS_v','lat','lon','time','stations')