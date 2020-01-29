%Este script ejecuta bajo entorno Matlab los programas
%de análisis espacial objetivo realizado por Damiá
%Gomis y Simón Ruiz del Institut Mediterrani d'Estudis
%AvanÇats (Barcelona).
%
%Entorno Matlab por Eugenio Fraile
%efraile@becarios.ulpgc.es
%
%15-10-02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%cargo los datos a procesar%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%El estudio se centrará en los perfiles de boyas 
%entre los meses: start y ending
start=4;%ej. 4=abril
ending=12;%ej. 9=septiembre

%cargo los datos de las boyas las 19 boyas
[temp,sal,sgth,pres,lat,lon,date]=...
   load_float_gyros_season(start,ending);

%inpath, siempre debe ser este (existe una organización
%de carpetas prefijada por los de IMEDEA
%todas las carpetas deben estar dentro de este path
inpath='D:\misdocumentos\Eugenio\Gyroscope\gyros0302\analisis_objetivo2\';

%ruta datos de salida (igual que el inpath)
outpath=inpath;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parámetros requeridos para la ejecución de los
%programas pero que no suelen cambiarse frecuentemente
%(al menos, para el análisis de los datos de una misma 
%campaña).

%Coordenadas de la esquina inferior izquierda
xlon1=-42;
xlat1=24;

%inclinación de la malla (en grados), deberá ser cero
%para una malla a lo largo de meridianos y paralelos
alfo=0;

%número de niveles verticales. Número de filas y
%columnas del nuevo grid.
NL=100;
NC=24;%como xarm=1, llegará hasta una lon=-19
NF=7;%como yarm=1, llegará hasta una lat=30

%tamaño de los brazos de malla horizontal (x,y): en
%grados, si la malla es a lo largo de meridianos y 
%y paralelos y en kilómetros si está inclinada
%debe existir un compromiso: 1/4*dm < Ax < 1/2*dm
%siendo:
% dm ==> distancia promedio entre estaciones (grados
%        si alfo=0, o si no, en km)
% Ax ==> el brazo de malla
xarm=1;
yarm=1;

%profundidad del primer nivel , intervalo vertical
%entre niveles y profundidad total.
zo=5;
zint=20;
ztot=2000;

%prefijo de los perfiles (s4) nombre estudio (s8)
suf='gyro';
studa='gyroscop';

%ejecutamos los programas
make_grid_data...
   (outpath,xlon1,xlat1,alfo,NL,NC,NF,xarm,yarm,zo,zint,suf,studa);

make_ctd_list(inpath,start,ending,outpath,suf,...
   temp,sal,sgth,pres,lat,lon,date);

%se ejecuta para la variable temperatura
run_analisis_obj(outpath);

make_station_grid(outpath,suf,zo,zint,ztot,NL,NC,NF);

%y ahora para variable salinidad
run_analisis_obj(outpath);

make_station_grid(outpath,suf,zo,zint,ztot,NL,NC,NF);

%finalmente se hace el postprocesado de los datos
%quedando en la carpeta mat los perfiles individuales
%de la periferia de la caja y uno global.
post_analisis_obj;
