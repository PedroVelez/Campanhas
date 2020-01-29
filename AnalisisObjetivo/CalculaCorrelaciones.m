warning off;close all;clear all

file='C:\Fisica\Proyectos\Tunibal\Tunibal2005\Hidrografia\CTD\Datos\Secciones\Total_Tbal2005';load(file)
data=sals;
pre=20;

indice=find(pres == pre);
zobs=data(indice,:)';
yobs=lats';
xobs=lons';
[zi]=polyfitn(1,1,xobs,yobs,zobs);
zobs=zobs-zi; %Solamente trabajo con anomalias

dint=1; %Intervalo de promedio en distancias
distmax=250; %Lag maximo sobre el que se realiza el calculo de correlaciones.
disti=3;
distf=250;

[corisom,coriso,dist] = croscor(xobs,yobs,zobs,dint,distmax);


%Quito valores cercano a lag=0 y para distancias muy largas.
dist(dist>distf)=NaN;
coriso=coriso(~isnan(dist));
dist=dist(~isnan(dist));

dist(dist<disti)=NaN;
coriso=coriso(~isnan(dist));
dist=dist(~isnan(dist));

dist=dist(~isnan(coriso));
coriso=coriso(~isnan(coriso));

para=lsqcurvefit('FCorGaus',[0.1 15],dist,coriso);gamma=para(1);scl=para(2);

disp(sprintf('     >Gamma=%5.3f scl=%5.3f',gamma,scl))

plot((0:1:max(dist))',f_correlacionobservada(para,(0:1:max(dist))'),'b','linewidth',2)
title(sprintf('Multivariable normalized autocovarianze \n Gamma=%5.3f scl=%5.3f',gamma,scl),'interpreter','none')
xlabel('Distance (km)');ylabel('Scale Correlation')
legend('Autocovarianza multivariable normalizada','Autocovarianza normalizada','Gausiana ajustada')
