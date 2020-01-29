Limpia

SVOpciones.titulo='RaProCan1804';
ZrangeU=[-30:5:30];
ZisoU=[-20 -10 10 20];
SVOpciones.figuras=[4 7];
SVOpciones.colormap=BlueWhiteRedColormap;


%% Seccion Norte
load Ra1804LadcpNorte
TituloSeccion='Norte';
SVOpciones.tipo='Zonal';
SVOpciones.stations=nstats;

SVOpciones.batymetry='BatimetriaRa1804Norte.mat';
SVOpciones.xtitulo='Longitude';
SVOpciones.ytitulo='Depth (m)';
SVData.xvar=longsladcp;
SVData.yvar=-sw_dpth(presiladcp,mean(latisladcp));

SVOpciones.xrange=[-18.5 -12.46];
SVOpciones.yrange1=[-400  -10];
SVOpciones.yrange2=[-4000 -400];

SVOpciones.ztitulo=strcat(TituloSeccion,' U-Velocidad');
SVData.zvar=uiladcp*100; %pasi a cm/s
SVOpciones.zrange=ZrangeU;     %
SVOpciones.ziso=ZisoU;       
figure;[ax1,ax2,ax3]=SeccionVertical(SVData,SVOpciones);

SVOpciones.ztitulo=strcat(TituloSeccion,' V-Velocidad');
SVData.zvar=viladcp*100; %pasi a cm/s
SVOpciones.zrange=ZrangeU;     
SVOpciones.ziso=ZisoU;       
figure;[ax1,ax2,ax3]=SeccionVertical(SVData,SVOpciones);

SVOpciones.ztitulo=strcat(TituloSeccion,' V(p)-Velocidad');
SVData.zvar=vpiladcp*100; %pasi a cm/s
SVOpciones.zrange=ZrangeU;     
SVOpciones.ziso=ZisoU;       
figure;[ax1,ax2,ax3]=SeccionVertical(SVData,SVOpciones);

%% Seccion Norte
load Ra1804Ladcp2Norte
TituloSeccion='Norte2';
SVOpciones.tipo='Zonal';
SVOpciones.stations=nstats;

SVOpciones.batymetry='BatimetriaRa1804Norte.mat';
SVOpciones.xtitulo='Longitude';
SVOpciones.ytitulo='Depth (m)';
SVData.xvar=longsladcp;
SVData.yvar=-sw_dpth(presiladcp,mean(latisladcp));

SVOpciones.xrange=[-18.5 -12.46];
SVOpciones.yrange1=[-400  -10];
SVOpciones.yrange2=[-4000 -400];

SVOpciones.ztitulo=strcat(TituloSeccion,' U-Velocidad');
SVData.zvar=uiladcp*100; %pasi a cm/s
SVOpciones.zrange=ZrangeU;     %
SVOpciones.ziso=ZisoU;       
figure;[ax1,ax2,ax3]=SeccionVertical(SVData,SVOpciones);

SVOpciones.ztitulo=strcat(TituloSeccion,' V-Velocidad');
SVData.zvar=viladcp*100; %pasi a cm/s
SVOpciones.zrange=ZrangeU;     
SVOpciones.ziso=ZisoU;       
figure;[ax1,ax2,ax3]=SeccionVertical(SVData,SVOpciones);

SVOpciones.ztitulo=strcat(TituloSeccion,' V(p)-Velocidad');
SVData.zvar=vpiladcp*100; %pasi a cm/s
SVOpciones.zrange=ZrangeU;     
SVOpciones.ziso=ZisoU;       
figure;[ax1,ax2,ax3]=SeccionVertical(SVData,SVOpciones);


%% Seccion Norte
load Ra1804LadcpLaGraciosa
TituloSeccion='LaGraciosa';
SVOpciones.tipo='m';
SVOpciones.stations=nstats;

SVOpciones.batymetry='BatimetriaRa1804LaGraciosa.mat';
SVOpciones.xtitulo='Latitude';
SVOpciones.ytitulo='Depth (m)';
SVData.xvar=latisladcp;
SVData.yvar=-sw_dpth(presiladcp,mean(latisladcp));

SVOpciones.xrange=[29.44 31.1];
SVOpciones.yrange1=[-400  -10];
SVOpciones.yrange2=[-3000 -400];

SVOpciones.ztitulo=strcat(TituloSeccion,' U-Velocidad');
SVData.zvar=uiladcp*100; %pasi a cm/s
SVOpciones.zrange=ZrangeU;     %
SVOpciones.ziso=ZisoU;       
figure;SeccionVertical(SVData,SVOpciones);

SVOpciones.ztitulo=strcat(TituloSeccion,' V-Velocidad');
SVData.zvar=viladcp*100; %pasi a cm/s
SVOpciones.zrange=ZrangeU;     
SVOpciones.ziso=ZisoU;       
figure;SeccionVertical(SVData,SVOpciones);

SVOpciones.ztitulo=strcat(TituloSeccion,' V(p)-Velocidad');
SVData.zvar=vpiladcp*100; %pasi a cm/s
SVOpciones.zrange=ZrangeU;     
SVOpciones.ziso=ZisoU;       
figure;[ax1,ax2,ax3]=SeccionVertical(SVData,SVOpciones);


%% CaboGhir
load Ra1804LadcpCaboGhir
TituloSeccion='CaboGhir';
SVOpciones.tipo='Zonal';
SVOpciones.stations=nstats;

SVOpciones.batymetry='BatimetriaRa1804CaboGhir.mat';
SVOpciones.xtitulo='Longitude';
SVOpciones.ytitulo='Depth (m)';
SVData.xvar=longsladcp;
SVData.yvar=-sw_dpth(presiladcp,mean(latisladcp));

SVOpciones.xrange=[-13.54 -10.0];
SVOpciones.yrange1=[-400  -10];
SVOpciones.yrange2=[-3000 -400];

SVOpciones.ztitulo=strcat(TituloSeccion,' U-Velocidad');
SVData.zvar=uiladcp*100; %pasi a cm/s
SVOpciones.zrange=ZrangeU;     %
SVOpciones.ziso=ZisoU;       
figure;[ax1,ax2,ax3]=SeccionVertical(SVData,SVOpciones);

SVOpciones.ztitulo=strcat(TituloSeccion,' V-Velocidad');
SVData.zvar=viladcp*100; %pasi a cm/s
SVOpciones.zrange=ZrangeU;     
SVOpciones.ziso=ZisoU;       
figure;[ax1,ax2,ax3]=SeccionVertical(SVData,SVOpciones);