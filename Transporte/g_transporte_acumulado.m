close all;clear all

figuras=[4 7];
iCapaSuperior=1:3;
iCapaIntermedia=4:6;
iCapaProfunda=7:12;

load ../../DatosCampanha
inpath='/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/2016_Raprocan1603/CTD/Datos/Mat/';

%stations=[1:51];
stations=[1:49];

%Lanzarote
load trans_masa_Lanzarote
mass_trans_lanzarote=mass_trans;
% load ../VelocidadLADCP/refvel_Raprocan1603_Lanzarote
% refvel(isnan(refvel)==1)=0;
% mass_trans_lanzarote_ref=mass_trans+mass.*[ones(12,1)*refvel];
mass_trans_lanzarote_ref=mass_trans;


%Norte
load trans_masa_Norte
mass_trans_norte=mass_trans;
% load ../VelocidadLADCP/refvel_Raprocan1603_Norte
% refvel(isnan(refvel)==1)=0;
% mass_trans_norte_ref=mass_trans+mass.*[ones(12,1)*refvel];
mass_trans_norte_ref=mass_trans;

%Oeste
load trans_masa_Oeste
mass_trans_oeste=mass_trans;
% load ../VelocidadLADCP/refvel_Raprocan1603_Oeste
% refvel(isnan(refvel)==1)=0;
% mass_trans_oeste_ref=mass_trans+mass.*[ones(12,1)*refvel];
mass_trans_oeste_ref=mass_trans;

% %Sur
 load trans_masa_Sur
 mass_trans_sur=mass_trans;
% load ../VelocidadLADCP/refvel_raprocan1603_Sur
% refvel(isnan(refvel)==1)=0;
% mass_trans_sur_ref=mass_trans+mass.*[ones(12,1)*refvel];
mass_trans_sur_ref=mass_trans;
 
mass_trans=    [fliplr(mass_trans_lanzarote)     ones(12,1)*0 fliplr(mass_trans_norte)     -mass_trans_oeste      -mass_trans_sur];
mass_trans_ref=[fliplr(mass_trans_lanzarote_ref) ones(12,1)*0 fliplr(mass_trans_norte_ref) -mass_trans_oeste_ref  -mass_trans_sur_ref];


% Tengo que llamar al get_stas porque necesito las posciones de las
% estaciones
[lat,lon,s,t,p,pt,sgth,gm,dh,pe]=get_stas_Raprocan1603(inpath,stations,10);
% Aqui calcula la distancia entre las estaciones
[dist,phaseangle]=sw_dist(lat,lon,'km');
dist=dist';
distcum=[0; cumsum(dist)];

% Estas son las distancias en las que empiezan cada transecto
distm(1)=distcum(10); % Empieza Lanzarote
distm(2)=distcum(11); % Acaba Lanzarote
distm(3)=distcum(24); % Fin transecto norte
distm(4)=distcum(30); % Fin transecto Oeste

%% Transporte acumulado
mass_trans_sup=sum(mass_trans(iCapaSuperior,:));
mass_trans_int=sum(mass_trans(iCapaIntermedia,:));
mass_trans_prof=sum(mass_trans(iCapaProfunda,:));
mass_trans_sup_ref=sum(mass_trans_ref(iCapaSuperior,:));
mass_trans_int_ref=sum(mass_trans_ref(iCapaIntermedia,:));
mass_trans_prof_ref=sum(mass_trans_ref(iCapaProfunda,:));
%Calculo el transporte acumulado
masscum_sup=[0; cumsum(mass_trans_sup')];
masscum_int=[0; cumsum(mass_trans_int')];
masscum_prof=[0; cumsum(mass_trans_prof')];
masscum_tot=masscum_sup+masscum_int+masscum_prof;
masscum_sup_ref=[0; cumsum(mass_trans_sup_ref')];
masscum_int_ref=[0; cumsum(mass_trans_int_ref')];
masscum_prof_ref=[0; cumsum(mass_trans_prof_ref')];
masscum_tot_ref=masscum_sup_ref+masscum_int_ref+masscum_prof_ref;

save(strcat(campanha,'_masa'))

%% Figuras
figure(1)
e=plot(distcum(1:10),masscum_sup_ref(1:10)/1e9,'r','linewidth',2.5);hold on;grid on;
plot(distcum(1:10),masscum_sup(1:10)/1e9,'r:','linewidth',2.5);

plot(distcum(11:end),masscum_sup_ref(11:end)/1e9,'r','linewidth',2.5);hold on;grid on;
plot(distcum(11:end),masscum_sup(11:end)/1e9,'r:','linewidth',2.5);

f=plot(distcum(1:10),masscum_int_ref(1:10)/1e9,'g','linewidth',2.5);
plot(distcum(1:10),masscum_int(1:10)/1e9,'g:','linewidth',2.5);

plot(distcum(11:end),masscum_int_ref(11:end)/1e9,'g','linewidth',2.5);
plot(distcum(11:end),masscum_int(11:end)/1e9,'g:','linewidth',2.5);

j=plot(distcum(1:10),masscum_prof_ref(1:10)/1e9,'b','linewidth',2.5);
plot(distcum(1:10),masscum_prof(1:10)/1e9,'b:','linewidth',2.5);

plot(distcum(11:end),masscum_prof_ref(11:end)/1e9,'b','linewidth',2.5);
plot(distcum(11:end),masscum_prof(11:end)/1e9,'b:','linewidth',2.5);

%Esto es para que me pinta una linea cuando empieza cada seccion
for i=1:size(distm,2)
    plot([distm(i) distm(i)], [-15 10],'linewidth',2,'color',[0.75 0.75 0.75])
end
axis([0 max(distcum) -6 5]);

title(sprintf('Acumulated Mass Transport %s',campanha),'fontsize',12);
xlabel('Distance (km)','fontsize',12)
ylabel('Mass Transport (10^9 kg s^{-1})','fontsize',12)
%legend_handles=[e f j k];
legend_handles=[e f j];
legend(legend_handles,sprintf('Surface-%4.2f',cuts(iCapaSuperior(2))), ...
    sprintf('%4.2f-%4.2f',cuts(iCapaIntermedia(1)),cuts(iCapaIntermedia(2))), ...
    sprintf('%4.2f-bottom',cuts(iCapaProfunda(1))), ... %'Net', ...
    3)
add_scale('t',distcum,stations);%pone el numero secuencial

orient tall;orient landscape;
CreaFigura(1,strcat(mfilename,campanha),figuras)

figure(2)
e=plot(distcum(1:10),masscum_sup_ref(1:10)/1e9,'r','linewidth',2.5);hold on;grid on;
plot(distcum(1:10),masscum_sup(1:10)/1e9,'r:','linewidth',2.5);
plot(distcum(11:end),masscum_sup_ref(11:end)/1e9,'r','linewidth',2.5);hold on;grid on;
plot(distcum(11:end),masscum_sup(11:end)/1e9,'r:','linewidth',2.5);
%Esto es para que me pinta una linea cuando empieza cada seccion
for i=1:size(distm,2)
    plot([distm(i) distm(i)], [-15 10],'linewidth',2,'color',[0.75 0.75 0.75])
end
axis([0 max(distcum) -12 12]);

title(sprintf('Acumulated Mass Transport %s',campanha),'fontsize',12);
xlabel('Distance (km)','fontsize',12)
ylabel('Mass Transport (10^9 kg s^{-1})','fontsize',12)
add_scale('t',distcum,stations);%pone el numero secuencial

figure(3)
f=plot(distcum(1:10),masscum_int_ref(1:10)/1e9,'g','linewidth',2.5);hold on;grid on
plot(distcum(1:10),masscum_int(1:10)/1e9,'g:','linewidth',2.5);

plot(distcum(11:end),masscum_int_ref(11:end)/1e9,'g','linewidth',2.5);
plot(distcum(11:end),masscum_int(11:end)/1e9,'g:','linewidth',2.5);

%Esto es para que me pinta una linea cuando empieza cada seccion
for i=1:size(distm,2)
    plot([distm(i) distm(i)], [-15 10],'linewidth',2,'color',[0.75 0.75 0.75])
end
axis([0 max(distcum) -9 3]);
title(sprintf('Acumulated mass transport %4.2f-%4.2f, %s',cuts(iCapaIntermedia(1)),cuts(iCapaIntermedia(2)),campanha),'fontsize',12);xlabel('Distance (km)','fontsize',12)
ylabel('Mass Transport (10^9 kg s^{-1})','fontsize',12)
add_scale('t',distcum,stations);%pone el numero secuencial

figure(4)
j=plot(distcum(1:10),masscum_prof_ref(1:10)/1e9,'b','linewidth',2.5);hold on;grid on
plot(distcum(1:10),masscum_prof(1:10)/1e9,'b:','linewidth',2.5);
plot(distcum(11:end),masscum_prof_ref(11:end)/1e9,'b','linewidth',2.5);
plot(distcum(11:end),masscum_prof(11:end)/1e9,'b:','linewidth',2.5);
%Esto es para que me pinta una linea cuando empieza cada seccion
for i=1:size(distm,2)
    plot([distm(i) distm(i)], [-15 10],'linewidth',2,'color',[0.75 0.75 0.75])
end
axis([0 max(distcum) -9 3]);
title(sprintf('% Acumulated Mass Transport 4.2f-bottom, %s',cuts(iCapaProfunda(1)),campanha),'fontsize',12);
xlabel('Distance (km)','fontsize',12)
ylabel('Mass Transport (10^9 kg s^{-1})','fontsize',12)
add_scale('t',distcum,stations);%pone el numero secuencial

