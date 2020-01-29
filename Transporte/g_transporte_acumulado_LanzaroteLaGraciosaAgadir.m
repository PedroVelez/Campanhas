close all;clear all

figuras=[4 7];

dd=10;
stations=[1:1:10 1101:1:1109 1016:-1:1001];

iCapaSuperior=1:3;
iCapaIntermedia=4:6;
iCapaProfunda=7:12;

load ../../DatosCampanha
inpath='/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/2017_Raprocan1710/CTD/Mat/';

%% Lanzarote
load trans_masa_Lanzarote
mass_trans_lanzarote=mass_trans;
%Anado el transporte debido a la velocidad en la capa de referencia.
load ../VelocidadLADCP/refvel_Raprocan1710_Lanzarote
refvel(isnan(refvel)==1)=0;
mass_trans_lanzarote_ref=mass_trans+mass.*[ones(12,1)*refvel];

%% La Graciosa
load trans_masa_LaGraciosa
mass_trans_LaGraciosa=mass_trans;
%Anado el transporte debido a la velocidad en la capa de referencia.
load ../VelocidadLADCP/refvel_Raprocan1710_LaGraciosa
refvel(isnan(refvel)==1)=0;
mass_trans_LaGraciosa_ref=mass_trans+mass.*[ones(12,1)*refvel];

%% Cabo Ghir
load trans_masa_CaboGhir
mass_trans_CaboGhir=mass_trans;
%Anado el transporte debido a la velocidad en la capa de referencia.
load ../VelocidadLADCP/refvel_Raprocan1710_CaboGhir
refvel(isnan(refvel)==1)=0;
mass_trans_CaboGhir_ref=mass_trans+mass.*[ones(12,1)*refvel];


%% Acumulo
mass_trans=    [-fliplr(mass_trans_lanzarote)     ones(12,1).*0 -fliplr(mass_trans_LaGraciosa)  mass_trans_CaboGhir];
mass_trans_ref=[-fliplr(mass_trans_lanzarote_ref) ones(12,1).*0 -fliplr(mass_trans_LaGraciosa_ref) mass_trans_CaboGhir_ref];

% Tengo que llamar al get_stas poruqe solo quiero una presion
[lat,lon,s,t,p,pt,sgth,gm,dh,pe]=get_stas_Raprocan1710(inpath,stations,dd);

% Aqui calcula la distancia entre las estaciones
[dist,phaseangle]=sw_dist(lat,lon,'km');
dist=dist';
distcum=[0; cumsum(dist)];

% Quito un
distcum=[distcum(1:10)' distcum(11:end)'-65]';

%Transporte
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

k=plot(distcum(1:10),masscum_tot_ref(1:10)/1e9,'k','linewidth',2.5);
plot(distcum(1:10),masscum_tot(1:10)/1e9,'k:','linewidth',2.5);

plot(distcum(11:end),masscum_tot_ref(11:end)/1e9,'k','linewidth',2.5);
plot(distcum(11:end),masscum_tot(11:end)/1e9,'k:','linewidth',2.5);

axis([0 max(distcum) -8 6]);

title(sprintf('Acumulated Mass Transport LanzaroteLaGraciosaCaboGhir %s',campanha),'fontsize',12);
xlabel('Distance (km)','fontsize',14)
ylabel('Mass Transport (10^9 kg s^{-1})','fontsize',14)
legend_handles=[e f j k];
legend(legend_handles,sprintf('Surface-%4.2f',cuts(iCapaSuperior(2))), ...
    sprintf('%4.2f-%4.2f',cuts(iCapaIntermedia(1)),cuts(iCapaIntermedia(2))), ...
    sprintf('%4.2f-bottom',cuts(iCapaProfunda(1))),'Net', ...
    'best')
add_scale('t',distcum,stations);%pone el numero secuencial

orient tall;orient landscape
CreaFigura(1,strcat(mfilename,campanha),figuras)