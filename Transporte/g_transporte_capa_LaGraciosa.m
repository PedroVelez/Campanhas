Limpia
%Grafico el transporte por capas y transecto

DC=load('../DatosCampanha');

load trans_masa_LaGraciosa
mass_trans_LaGraciosa=mass_trans;
%load ../VelocidadLADCP/refvel_Raprocan1610_Norte
%refvel(isnan(refvel)==1)=0;
%mass_trans_norte_ref=mass_trans+mass.*[ones(12,1)*refvel];
mass_trans_norte_ref=mass_trans;

ha=plot(sum(mass_trans_norte_ref')/1e9,(1:12),'bo-','markersize',4,'linewidth',2.25);hold on
hb=plot(sum(mass_trans_LaGraciosa')    /1e9,(1:12),'bo:','markersize',4,'linewidth',2.25);hold on
plot([0 0],[0 15],'r--','linewidth',1);grid on;set(gca,'YDir','reverse')

YTickLabel=['Surfac';'26.440';'26.850';'27.162';'27.380';'27.620';'27.820';'27.922';'27.975';'28.008';'28.044';'28.072';'28.099';'28.110';'bottom'];
set(gca,'YTick',[0:1:14]+.5,'YTickLabel',YTickLabel) ;
xlabel('Sv','fontsize',14)
ylabel('\gamma_n (kg m^{-3})','fontsize',14)
title('Initial Mass Transport La Graciosa','fontsize',14)
legend_handles=[ha hb];
legend(legend_handles,'Northern Transect referenced','Northern Transect')
axis([-3 3 0.5 14.5])

CreaFigura(gcf,strcat(mfilename,DC.campanha),[4 7])
