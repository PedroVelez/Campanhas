% MapaTodasEstaciones
%		Dibuja el mapa con las posiciones de todas las estaciones
%		que hay en el subdirectorio ./mat
%
%		v1.0 25 Junio 2003 - Instituto Espa?ol de Oceanografia
%		v1.1 25 Junio 2004 - Instituto Espa?ol de Oceanografia
%		v1.2 25 Junio 2005 - Instituto Espa?ol de Oceanografia
%			Hacemos cambios para que no de errores cuando no hay ficheros de batimetria y/o costa
%
Limpia

Config.Batimetria=1;        %[1/0] for batythemtry
Config.BatimetriaIsobaths=[-1000 -2000 -4000];%Isobaths to contour
Config.BatimetriaIsobathsLabel=[1 1 1];%Isobaths to label
Config.BatimetriaColor=0;   %[1/0] for backgroung batymetry.

Config.StationTicks=5;    %Cada cuanto se etiquetan las estaciones

Config.StationSpecialMarks1=[];
Config.StationSpecialMarks1Color='r';
Config.StationSpecialMarks1Ticks=0;

Config.OutputFigures=[4 7];

Config.ImagenSatelite=0;     %Si hay imagen por satelite la incluye de fondo. [1]MADT [2]SST [3] Manual
Config.ImagenSateliteTitulo='Chl 9 Octubre';
% Config.ImagenSateliteFile='201310151405MSN19';

%% Begin
GlobalDS=load(fullfile('..','DatosCampanha'));
if isfield(GlobalDS,'campanha')==0; GlobalDS.campanha='';end
if isfield(GlobalDS,'filebat')==0; GlobalDS.filebat='';end
if isfield(GlobalDS,'filecosta')==0; GlobalDS.filecosta='';end
if isfield(GlobalDS,'filemalla')==0; GlobalDS.filemalla='';end
if isfield(GlobalDS,'filemtopo')==0; GlobalDS.filemtopo='';end

m_proj('Mercator','long',[GlobalDS.lon_min GlobalDS.lon_max],'lat',[GlobalDS.lat_min GlobalDS.lat_max]);

if Config.ImagenSatelite==3
    DCHL=load('/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/Administracion/2017_Raprocan1710/PlanCampanha/Modis_chla_Raprocan1017.mat');
    [CHlon,CHlat]=meshgrid(DCHL.lon+360,DCHL.lat);
    L=squeeze(nanmean(DCHL.chla_all,1));
    [X,Y]=m_ll2xy(CHlon,CHlat,'clip','on');
    contourf(X,Y,L,.1:0.005:.3,'edgecolor','none');hold on
    caxis([.1 .30])
end

if Config.Batimetria==1
    if ~isempty(GlobalDS.filebat)
        BAT=load(GlobalDS.filebat);
        Config.BatimetriaIsobaths=sort(Config.BatimetriaIsobaths,'descend');
        cbiso=linspace(0.5,1,length(Config.BatimetriaIsobaths));
        for iiso=1:length(Config.BatimetriaIsobaths)
            [C,h]=m_contour(BAT.batylon,BAT.batylat,BAT.elevations,[Config.BatimetriaIsobaths(iiso) Config.BatimetriaIsobaths(iiso)],'color',cbiso(iiso)*[1 1 1]);hold on
            if Config.BatimetriaIsobathsLabel(iiso)==1
                clabel(C,h,'FontSize',9,'LabelSpacing',500,'Color',cbiso(iiso)*[1 1 1])
            end
        end
    end
end

if Config.BatimetriaColor==1 || Config.ImagenSatelite>0
    mfc=[0.85 0.85 0.85]; %Marker Face Color
    mec='k'; %Marker edge color
    ctrack=[0.85,0.85,0.85];
else
    mfc=[0.45 0.45 0.45];
    mec='k';
    ctrack=[0.25,0.25,0.25];
end

files=dir(strcat('../../CTD',filesep,'Mat',filesep,'*.mat'));
fechas=[];lon=[];lat=[];nstat=[];

%leo los datos
for ifile=1:size(files,1)
    file=strcat('../../CTD',filesep,'Mat',filesep,files(ifile).name);
    if ~isempty(regexp(file,strcat('.',filesep,'Mat',filesep,GlobalDS.campanhacode,'_(\w+).mat'),'ignorecase','once'))
        fprintf('%s \n',file)
        EstData=load(file);
        fechas=[fechas datenum(EstData.gtime)];
        lon=[EstData.long+360 lon];
        lat=[EstData.lati lat];
        nstat=[EstData.nstat nstat];
    end
end
%ordeno los datos por fecha
[fechas,IA]=sort(fechas);
lon=lon(IA);
lat=lat(IA);
nstat=nstat(IA);

for iest=1:length(nstat)
    m_line(lon(iest),lat(iest),'marker','o','markersize',6,'Markeredgecolor',mec,'MarkerFaceColor',mfc);
    if Config.StationTicks>=1
        if ceil(iest/Config.StationTicks)==iest/Config.StationTicks
            m_text(lon(iest),lat(iest)+0.2,num2str(nstat(iest)),'color',mec,'Fontsize',11,'HorizontalAlignment','center','VerticalAlignment','top','fontname','cambria');
            m_line(lon(iest),lat(iest),'marker','o','markersize',8,'Markeredgecolor',mec,'MarkerFaceColor',mfc);
        end
    end
end
if ~isempty(GlobalDS.filecosta)
    m_usercoast(GlobalDS.filecosta,'patch',[.7 .6 .4,],'edgecolor',[.7 .6 .4,]);
end
m_grid('linestyle','none','fontsize',12,'fontname','cambria')
title(sprintf('%s %d estaciones, [%s - %s].\n',strcat(GlobalDS.campanha),length(fechas),datestr(min(fechas),1),datestr(max(fechas),1)),'fontsize',12,'fontname','cambria');

CreaFigura(1,strcat(mfilename,GlobalDS.campanha),Config.OutputFigures)