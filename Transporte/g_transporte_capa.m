close all;clear all
%Grafico el transporte por capas y transecto

campanha='Raprocan1504';


load ../../DatosCampanha


%Lanzarote
load trans_masa_Lanzarote
mass_trans_Lanzarote=mass_trans;
% load ../VelocidadLADCP/refvel_Raprocan1510_Lanzarote
% refvel(isnan(refvel)==1)=0;
% mass_trans_Lanzarote_ref=mass_trans+mass.*[ones(12,1)*refvel];
mass_trans_Lanzarote_ref=mass_trans;

%Norte
load trans_masa_Norte
mass_trans_Norte=mass_trans;
% load ../VelocidadLADCP/refvel_Raprocan1510_Norte
% refvel(isnan(refvel)==1)=0;
% mass_trans_Norte_ref=mass_trans+mass.*[ones(12,1)*refvel];
mass_trans_Norte_ref=mass_trans;

mass_trans_Norte=    [fliplr(mass_trans_Lanzarote)     ones(12,1)*0 fliplr(mass_trans_Norte)];
mass_trans_Norte_ref=[fliplr(mass_trans_Lanzarote_ref) ones(12,1)*0 fliplr(mass_trans_Norte_ref)];

%Oeste
load trans_masa_Oeste
mass_trans_Oeste=mass_trans;
% load ../VelocidadLADCP/refvel_Raprocan1510_Oeste
% refvel(isnan(refvel)==1)=0;
% mass_trans_Oeste_ref=mass_trans+mass.*[ones(12,1)*refvel];
mass_trans_Oeste_ref=mass_trans;

% Sur
load trans_masa_Sur
mass_trans_Sur=mass_trans;
% load ../VelocidadLADCP/refvel_raprocan1510_Sur
% refvel(isnan(refvel)==1)=0;
% mass_trans_Sur_ref=mass_trans+mass.*[ones(12,1)*refvel];
mass_trans_Sur_ref=mass_trans;

%Total (Sigue el convenio de la caja).
mass_trans_Total=    [mass_trans_Norte     -mass_trans_Oeste     -mass_trans_Sur];
mass_trans_Total_ref=[mass_trans_Norte_ref -mass_trans_Oeste_ref -mass_trans_Sur_ref];

figure
hN=plot(sum(mass_trans_Norte_ref')/1e9,(1:12),'ro-','markersize',4,'linewidth',2.25);hold on
plot(sum(mass_trans_Norte')    /1e9,(1:12),'ro:','markersize',4,'linewidth',2.25);hold on

hO=plot(sum(mass_trans_Oeste_ref')/1e9,(1:12),'go-','markersize',4,'linewidth',2.25);hold on
plot(sum(mass_trans_Oeste')    /1e9,(1:12),'go:','markersize',4,'linewidth',2.25);hold on

hS=plot(sum(mass_trans_Sur_ref')/1e9,(1:12),'bo-','markersize',4,'linewidth',2.25);hold on
plot(sum(mass_trans_Sur')    /1e9,(1:12),'bo:','markersize',4,'linewidth',2.25);hold on

hT=plot(sum(mass_trans_Total_ref')/1e9,(1:12),'ko-','markersize',4,'linewidth',2.25);hold on
plot(sum(mass_trans_Total')    /1e9,(1:12),'ko:','markersize',4,'linewidth',2.25);hold on

plot([0 0],[0 15],'r--','linewidth',2);grid on;set(gca,'YDir','reverse')

set(gca,'YTick',[0:1:13]+.5,'YTickLabel','Surface|26.440|26.850|27.162|27.380|27.620|27.820|27.922|27.975|28.008|28.044|28.072|28.099|bottom') ;
xlabel('Sv','fontsize',14)
ylabel('\gamma_n (kg m^{-3})','fontsize',14)
title('Initial Mass Transport','fontsize',14)
legend_handles=[hN hO hS hT];
legend(legend_handles,'North','West','South','Total')
axis([-3 3 0.5 13.5])

CreaFigura(gcf,strcat(mfilename,campanha),7)

figure
subplot(1,2,1)
hN=plot(sum(mass_trans_Norte')/1e9,(1:12),'ro-','markersize',4,'linewidth',2.25);hold on
hO=plot(sum(mass_trans_Oeste')/1e9,(1:12),'go-','markersize',4,'linewidth',2.25);hold on
hS=plot(sum(mass_trans_Sur')/1e9,(1:12),'bo-','markersize',4,'linewidth',2.25);hold on
hT=plot(sum(mass_trans_Total')/1e9,(1:12),'ko-','markersize',4,'linewidth',2.25);hold on
plot([0 0],[0 15],'r--','linewidth',2);grid on;set(gca,'YDir','reverse')
set(gca,'YTick',[0:1:13]+.5,'YTickLabel','Surface|26.440|26.850|27.162|27.380|27.620|27.820|27.922|27.975|28.008|28.044|28.072|28.099|bottom') ;
xlabel('Sv','fontsize',14)
ylabel('\gamma_n (kg m^{-3})','fontsize',14)
title('Initial Mass Transport','fontsize',14)
legend_handles=[hN hO hS hT];
legend(legend_handles,'North','West','South','Total')
axis([-3 3 0.5 13.5])

subplot(1,2,2)
hN=plot(sum(mass_trans_Norte_ref')/1e9,(1:12),'ro-','markersize',4,'linewidth',2.25);hold on
hO=plot(sum(mass_trans_Oeste_ref')/1e9,(1:12),'go-','markersize',4,'linewidth',2.25);hold on
hS=plot(sum(mass_trans_Sur_ref')/1e9,(1:12),'bo-','markersize',4,'linewidth',2.25);hold on
hT=plot(sum(mass_trans_Total_ref')/1e9,(1:12),'ko-','markersize',4,'linewidth',2.25);hold on
plot([0 0],[0 15],'r--','linewidth',2);grid on;set(gca,'YDir','reverse')
set(gca,'YTick',[0:1:13]+.5,'YTickLabel','Surface|26.440|26.850|27.162|27.380|27.620|27.820|27.922|27.975|28.008|28.044|28.072|28.099|bottom') ;
xlabel('Sv','fontsize',14)
ylabel('\gamma_n (kg m^{-3})','fontsize',14)
title('Initial Mass Transport LADCP','fontsize',14)
legend_handles=[hN hO hS hT];
legend(legend_handles,'North','West','South','Total')
axis([-3 3 0.5 13.5])

