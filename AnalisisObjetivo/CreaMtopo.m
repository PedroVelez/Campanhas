%Este programa crea una matriz con los valores
%interpolados de la batimetria en la malla del analisis objetivo
%Esta matriz se usara despues como mascara.
%v 1.0 - 25 Junio 2004 IEO
close all;clear all

%Area de analisis
%Area geografica
lat_min=35; lat_max=37.5;
lon_min=-6; lon_max=0;

%Area de analisis
xnew= 0.5:0.050: 5.05;%Mallado en x
ynew=37.8:0.050:40.40;%Mallado en y

filebat='Tunibalbat.mat'; %Ficheo con la batimetria
filemtopo='MtopoTbal0050'; %Fichero de salida de la matriz H

%-------------------------
[X,Y]=meshgrid(xnew,ynew);
load(filebat);
for i=1:size(X,1)
   for j=1:size(X,2)
      H(i,j)=elevations(locate(vlat,Y(i,j)),locate(vlon,X(i,j)));
   end
end   
keep H filemtopo
save(filemtopo,'H')
