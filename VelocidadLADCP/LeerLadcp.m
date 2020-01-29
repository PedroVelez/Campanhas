Limpia

DataDir='/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/2018_Raprocan1804/LADCP/Visbeck/profiles/';

DC=load('../DatosCampanha');
load(DC.filebat)

stations=[1:1:24 901:914 1101:1:1110 101:108];
lev=8:8:1664;
iest=0;
fprintf('Reading data \n  ')
for i2=1:length(stations)
    station=stations(i2);
    FileData=sprintf('%sra1804_%03d.mat',DataDir,station);
    if exist(FileData,'file')> 0
        fprintf('%2.2i, ',station)
        iest=iest+1;
        load(FileData)
        lons(iest)=dr.lon;
        lats(iest)=dr.lat;
        Data(iest).u=dr.u;
        Data(iest).v=dr.v;
        Data(iest).z=dr.z;
        Fechas(iest)=datenum(dr.date(1),dr.date(2),dr.date(3),dr.date(4),dr.date(5),dr.date(6));
        nstats(iest)=station;
    else
        fprintf('>>>> Falta el %sra1710_%03d.mat\n',DataDir,station);
    end
    
end

%% Agrupo a niveles homogeneos
fprintf('\n Averaging \n  ',station)
for i2=1:length(nstats)
    up=Data(i2).u;
    vp=Data(i2).v;
    zp=Data(i2).z;
    ulev(1:length(lev),i2)=interp1(zp,up,lev);
    vlev(1:length(lev),i2)=interp1(zp,vp,lev);
end

%% Figuras
EsVe=3;
xlege=-18;
ylege=31;
ulege=0.15;


figure
quiver(lev.*0,-lev,EsVe*nanmean(ulev'),EsVe*nanmean(vlev'),0)
title('Mean Velocity profile')
ylabel('Depth')
xlabel('cm/s')

figure
z=[10 88];
m_proj('Mercator','long',[DC.lon_min DC.lon_max],'lat',[DC.lat_min DC.lat_max]);
zi=Locate(lev,z(1));
zf=Locate(lev,z(2));
m_contour(batylon,batylat,elevations,[-4000 -4000],'color',[0.45 0.45 0.45]);hold on
m_contour(batylon,batylat,elevations,[-2000 -2000],'color',[0.55 0.55 0.55]);
m_contour(batylon,batylat,elevations,[-1000 -1000],'color',[0.65 0.65 0.65]);
m_contour(batylon,batylat,elevations,[ -500  -500],'color',[0.75 0.75 0.75]);
m_contour(batylon,batylat,elevations,[ -250  -250],'color',[0.85 0.85 0.85]);
m_usercoast(DC.filecosta,'patch',[.7 .6 .4,],'edgecolor',[.7 .6 .4,]);
title(sprintf('LADCP en %s between %03d-%03d m',DC.campanha,lev(zi),lev(zf)))
m_grid('box','on','ticklength',0.02,'linestyle','none')
X=[xlege lons]+360;
Y=[ylege lats];
U=EsVe*[ulege nanmean(ulev(zi:zf,:))];
V=EsVe*[ulege nanmean(vlev(zi:zf,:))];
m_quiver(X,Y,U,V,0,'k');hold on
m_text(xlege+360,ylege,strcat(num2str(ulege*100),' cm/s'),'color','k')
CreaFigura(gcf,sprintf('%s_%s_%04d_%04d',mfilename,DC.campanha,lev(zi),lev(zf)),[7 4]);


figure
z=[300 540];
m_proj('Mercator','long',[DC.lon_min DC.lon_max],'lat',[DC.lat_min DC.lat_max]);
zi=Locate(lev,z(1));
zf=Locate(lev,z(2));
m_contour(batylon,batylat,elevations,[-4000 -4000],'color',[0.45 0.45 0.45]);hold on
m_contour(batylon,batylat,elevations,[-2000 -2000],'color',[0.55 0.55 0.55]);
m_contour(batylon,batylat,elevations,[-1000 -1000],'color',[0.65 0.65 0.65]);
m_contour(batylon,batylat,elevations,[ -500  -500],'color',[0.75 0.75 0.75]);
m_contour(batylon,batylat,elevations,[ -250  -250],'color',[0.85 0.85 0.85]);
m_usercoast(DC.filecosta,'patch',[.7 .6 .4,],'edgecolor',[.7 .6 .4,]);
title(sprintf('LADCP en %s between %03d-%03d m',DC.campanha,lev(zi),lev(zf)))
m_grid('box','on','ticklength',0.02,'linestyle','none')
X=[xlege lons]+360;
Y=[ylege lats];
U=EsVe*[ulege nanmean(ulev(zi:zf,:))];
V=EsVe*[ulege nanmean(vlev(zi:zf,:))];
m_quiver(X,Y,U,V,0,'k');hold on
m_text(xlege+360,ylege,strcat(num2str(ulege*100),' cm/s'),'color','k')
CreaFigura(gcf,sprintf('%s_%s_%04d_%04d',mfilename,DC.campanha,lev(zi),lev(zf)),[7 4]);


figure
z=[900 1200];
m_proj('Mercator','long',[DC.lon_min DC.lon_max],'lat',[DC.lat_min DC.lat_max]);
zi=Locate(lev,z(1));
zf=Locate(lev,z(2));
m_contour(batylon,batylat,elevations,[-4000 -4000],'color',[0.45 0.45 0.45]);hold on
m_contour(batylon,batylat,elevations,[-2000 -2000],'color',[0.55 0.55 0.55]);
m_contour(batylon,batylat,elevations,[-1000 -1000],'color',[0.65 0.65 0.65]);
m_contour(batylon,batylat,elevations,[ -500  -500],'color',[0.75 0.75 0.75]);
m_contour(batylon,batylat,elevations,[ -250  -250],'color',[0.85 0.85 0.85]);
m_usercoast(DC.filecosta,'patch',[.7 .6 .4,],'edgecolor',[.7 .6 .4,]);
title(sprintf('LADCP en %s between %03d-%03d m',DC.campanha,lev(zi),lev(zf)))
m_grid('box','on','ticklength',0.02,'linestyle','none')
X=[xlege lons]+360;
Y=[ylege lats];
U=EsVe*[ulege nanmean(ulev(zi:zf,:))];
V=EsVe*[ulege nanmean(vlev(zi:zf,:))];
m_quiver(X,Y,U,V,0,'k');hold on
m_text(xlege+360,ylege,strcat(num2str(ulege*100),' cm/s'),'color','k')
CreaFigura(gcf,sprintf('%s_%s_%04d_%04d',mfilename,DC.campanha,lev(zi),lev(zf)),[7 4]);
