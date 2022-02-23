Limpia

ZisoTem=[5 10 15 20];
ZisoSal=[36.4 36.2 35.8 35.4 35.2];
ZisoGam=[26.44 26.85 27.162 27.38 27.62 27.82 27.922 27.975 28.008 28.044 28.072 28.0986 28.11 28.1295];
SVOpciones.figuras=[4 7];

load palette
%%
DC=load('../DatosCampanha');
SVOpciones.filecosta=DC.filecosta;
SVOpciones.lon_min=DC.lon_min;
SVOpciones.lon_max=DC.lon_max;
SVOpciones.lat_min=DC.lat_min;
SVOpciones.lat_max=DC.lat_max;
SVOpciones.titulo=DC.campanha;

%% Seccion Norte
if 1
    Data=load('Ra2012Norte');
    IndEst=1:length(Data.longs);
    SVOpciones.TituloSeccion='Norte';
    SVOpciones.tipo='Zonal';
    SVOpciones.stations=Data.nstats(IndEst);
    
    SVData.lon=Data.longs(:,IndEst);
    SVData.lat=Data.latis(:,IndEst);
    
    SVOpciones.batymetry='BatimetriaRa1903Norte.mat';
    SVOpciones.xtitulo='Longitude';
    SVOpciones.ytitulo='Depth (m)';
    SVData.xvar=Data.longs(:,IndEst);
    SVData.yvar=-sw_dpth(Data.press,mean(Data.latis(:,IndEst)));
    
    SVOpciones.xrange=[-18.5 -12.46];
    SVOpciones.yrange1=[-400  -10];
    SVOpciones.yrange2=[-4000 -400];
    
    figure;SeccionVerticalMapa(SVData,SVOpciones);
    
    SVOpciones.ztitulo=strcat(SVOpciones.TituloSeccion,' Temperatura');
    SVData.zvar=sw_ptmp(Data.salts(:,IndEst),Data.temps(:,IndEst),Data.press,0);
    SVOpciones.zrange=0;
    SVOpciones.ziso=ZisoTem;
    SVOpciones.colormap=WocePtmpColormap;
    %load palette
    SVOpciones.colormap=cmocean('thermal');
    figure;SeccionVertical(SVData,SVOpciones);
    
end
SVOpciones.ztitulo=strcat(SVOpciones.TituloSeccion,' Salinidad');
SVData.zvar=Data.salts(:,IndEst);
SVOpciones.ziso=ZisoSal;
%SVOpciones.colormap=WoceSaltColormap;
SVOpciones.colormap=cmocean('haline');
figure;SeccionVertical(SVData,SVOpciones);

SVOpciones.ztitulo=strcat(SVOpciones.TituloSeccion,' Gamma');
SVData.zvar=Data.gamas(:,IndEst);
SVOpciones.zrange=0;
SVOpciones.ziso=ZisoGam;
SVOpciones.colormap=WoceGammanColormap;
SVOpciones.colormap=palette;
figure;SeccionVertical(SVData,SVOpciones);

% %% Seccion LP
% Data=load('Ra2012LP');
% IndEst=1:length(Data.longs);
% SVOpciones.TituloSeccion='Lanzarote Passage';
% SVOpciones.tipo='Zonal';
% SVOpciones.stations=Data.nstats(IndEst);
%
% SVOpciones.batymetry='BatimetriaRa1903LP.mat';
% SVOpciones.xtitulo='Longitude';
% SVOpciones.ytitulo='Depth (m)';
% SVData.xvar=Data.longs(:,IndEst);
% SVData.yvar=-sw_dpth(Data.press,mean(Data.latis(:,IndEst)));
% SVData.lon=Data.longs(:,IndEst);
% SVData.lat=Data.latis(:,IndEst);
%
% SVOpciones.xrange=[-13.85 -12.47];
% SVOpciones.yrange1=[-400  -10];
% SVOpciones.yrange2=[-1400 -400];
%
% SVOpciones.ztitulo=strcat(SVOpciones.TituloSeccion,' Temperatura');
% SVData.zvar=sw_ptmp(Data.salts(:,IndEst),Data.temps(:,IndEst),Data.press,0);
% SVOpciones.zrange=0;
% SVOpciones.ziso=ZisoTem;
% SVOpciones.colormap=WocePtmpColormap;
% SVOpciones.colormap=palette;
% figure;SeccionVertical(SVData,SVOpciones);
%
% SVOpciones.ztitulo=strcat(SVOpciones.TituloSeccion,' Salinidad');
% SVData.zvar=Data.salts(:,IndEst);
% SVOpciones.ziso=ZisoSal;
% SVOpciones.colormap=WoceSaltColormap;
% figure;SeccionVertical(SVData,SVOpciones);
%
% SVOpciones.ztitulo=strcat(SVOpciones.TituloSeccion,' Gamma');
% SVData.zvar=Data.gamas(:,IndEst);
% SVOpciones.zrange=0;
% SVOpciones.ziso=ZisoGam;
% SVOpciones.colormap=WoceGammanColormap;
% figure;SeccionVertical(SVData,SVOpciones);
%
%

