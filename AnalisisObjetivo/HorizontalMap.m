clc;close all;clear all

load('..\..\..\DatosCampanha')
file='..\..\datos\Secciones\Total';
file='';

if strcmp(file,'')
    [file,pathfile]=uigetfile ( '..\..\datos\Secciones\*.mat', 'Escoge el ficheros de datos' );
    file=strcat(pathfile,file(1:findstr(file,'.mat')-1));
end    
load(strcat(file,'.mat')) ;%No poner .mat

profundidad=-25;%Profundidad donde se escoge el analisis. Si <0 ask for depth
automatico=1; 	%Variable para hacer contornos automaticos o no
etiquetas=2;  	%Variable para poner etiquetas en algunos contornos

disp('>>>>> (1) Salinity')
disp('>>>>> (2) Temperature')
disp('>>>>> (3) Density')
disp('>>>>> (4) Dynamic heigh and geostrophic velocity')
disp('>>>>> (5) Salinity and geostrophic velocity')
disp('>>>>> (6) Oxigen')

ivariable=input('>>>>>Which variable?');
switch ivariable
    case   1
        variable='Salinity';
        data=sals;        
        ies=0.05;   %intevalo de contorno
        escala=36.8:ies:38.0;	%Valores de contorno
        labelsetiquetas=[36.5,37,37.5,38,38.5]; %Etiquetas a poner
    case 2
        variable='Temperature';data=tems;
        ies=0.1;		%Intevalo de contorno
        escala=26.5:ies:28.5;	%Valores de contorno
        labelsetiquetas=[19,19.5,20,20.5]; %Etiquetas a poner
    case 3
        variable='Potential Density';data=pdes;
        ies=0.1;			%Intevalo de contorno
        escala=26.5:ies:28.5;	%Valores de contorno
        labelsetiquetas=[27,27.5,28,28.5]; %Etiquetas a poner
    case 4
        load(strcat(file,'Dyn600')); 
        variable='Geopotencial';data=10.*dyns; %paso el geopotencial a cmDyns
        fvector=1;		%Variable para poner vectores o no
        ies=0.5; 	%intevalo de contorno
        escala=-42.5:ies:-23.0; %Valores de contorno
        labelsetiquetas=[-2,0,2]; %Etiquetas a poner
    case 5
        variable='SalinityVec';data=sals;
        load(strcat(file,'Dyn600')); 
        fvector=1;		%Variable para poner vectores o no
        ies=0.05; 		%intevalo de contorno
        escala=36.7:ies:38.2;	%Valores de contorno
        labelsetiquetas=[36.5,37,37.5,38,38.5]; %Etiquetas a poner
    case 6
        variable='Oxygen';
        sato2=sw_sato2(sals,tems);
        data=oxys./sato2;
    otherwise
        disp('>>>>> Wrong variable');return
end

figuras=4;    	%Variable para hacer figura en 1 .jpeg, 2 .emf,3 .ps y 4.png
cl=jet;			%Escala de color
marcas=1;		%Variable para poner puntos con el valor real
msize=5; 		%Tamaño de los puntos
bati=1; 		%Variable para añadir topografia
batiprof=[-1000 -100];		%Contorno de batimetria.
mascara=1;		%Variable para usar las mascara
titlef=1;
subtitlef=1;	%Pone un subtitulo con los campos maximos/minimos observados

%Parametros del analisis
oberrormax=.3;			%Valor del analysis noise-to-signal ratio a partir del cual ya no considero el analisis
eta=0.05;
corlen=0.2;%.3 			%Escala de Correlación
merror=1;			%Variable para poner o no el mapa de errores

%Parametros para velocidades
munan=3;        %Velocidades menores que munan no se representan
dtrend=0;   	%Variable para quitar la media al campo observado de Dyn
fac=0.006;    	%Factor de escala para los vectores
ns=1; 			%Number of velocityvectors to skip

%Llamo a la subroutina de analisis y contorno
AnalisisYContornoHori