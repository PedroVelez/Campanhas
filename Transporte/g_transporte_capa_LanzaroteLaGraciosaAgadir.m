%Grafico el transporte por capas y transecto
Limpia

load ../../DatosCampanha

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


%% graficos
ha=plot(sum(mass_trans_ref')/1e9,(1:12),'bo-','markersize',4,'linewidth',2.25);hold on
hb=plot(sum(mass_trans')    /1e9,(1:12),'bo:','markersize',4,'linewidth',2.25);hold on
plot([0 0],[0 15],'r--','linewidth',1);grid on;set(gca,'YDir','reverse')

YTickLabel=['Surfac';'26.440';'26.850';'27.162';'27.380';'27.620';'27.820';'27.922';'27.975';'28.008';'28.044';'28.072';'28.099';'28.110';'bottom'];
set(gca,'YTick',[0:1:14]+.5,'YTickLabel',YTickLabel) ;

xlabel('Sv','fontsize',14)
ylabel('\gamma_n (kg m^{-3})','fontsize',14)
title('Initial Mass Transport LanzaroteLaGraciosaCaboGhir','fontsize',14)
legend_handles=[ha hb];
legend(legend_handles,'LADCP referenced','No-motion')
axis([-3 3 0.5 14.5])

CreaFigura(gcf,strcat(mfilename,campanha),[4 7])
