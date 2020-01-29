%DatosCampanha
%Fichero donde se almacena toda la informacion de la campanha en curso
%Usado por todos los programas del paquete HidroIEO
%
% 1.0 pedro.velez@ieo.es - 26 Junio 2004
close all;clear all;clc

campanha='Raprocan1903'; % Cruise name
campanhacode='Ra1903';   % Short cruise name

CruiseDir='/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/2019_Raprocan1903';
CruiseProyection='mercator';

filebat=strcat('CanaryIslandsBat.mat');  	% Fichero con la batimetria de la zona
filecosta=strcat('CanaryIslandsCoast');	% Fichero con la costa de la zona

%Geographical box
lat_min=24.25; lat_max=31.5;
lon_min=360-19; lon_max=360-8;

% Exteneded Geographical box
lat_min_Ext=22;   lat_max_Ext=34;
lon_min_Ext=360-24;  lon_max_Ext=360-8;

% Isolines to labels
BatimetriaIsobaths=[-1000 -2000 -4000 -5000];%Isobaths to contour
BatimetriaIsobathsLabel=[1 1 1 1];%Isobaths to label

%Place for titles
XL1=360-11.5;
YL1=25.75;
XL2=360-11.5;
YL2=24.75;

%Place for title Extended
XL1_Ext=360-12;
YL1_Ext=24;
XL2_Ext=360-9;
YL2_Ext=22.5;

%% figure
figure
m_proj(CruiseProyection,'lon',[lon_min_Ext lon_max_Ext],'lat',[lat_min_Ext lat_max_Ext])
BAT=load(filebat);
[CS,CH]=m_contourf(BAT.batylon,BAT.batylat,BAT.elevations, ...
    [-6000:100:0],'edgecolor','none');hold on

cbiso=linspace(0.5,1,length(BatimetriaIsobaths));
for iiso=1:length(BatimetriaIsobaths)
    [C,h]=m_contour(BAT.batylon,BAT.batylat,BAT.elevations, ...
        [BatimetriaIsobaths(iiso) BatimetriaIsobaths(iiso)], ...
        'color',cbiso(iiso)*[1 1 1]);hold on
    if BatimetriaIsobathsLabel(iiso)==1
        clabel(C,h,'FontSize',9,'LabelSpacing',500,'Color',cbiso(iiso)*[1 1 1])
    end
end
colormap(m_colmap('blues'));
m_usercoast(filecosta,'patch',[.7 .6 .4,],'edgecolor',[.7 .6 .4,]);hold on
m_line([lon_min lon_max lon_max lon_min lon_min],[lat_min lat_min lat_max lat_max lat_min], ...
    'color','k','linewidth',2)
m_grid('linestyle','none')
title(sprintf('%s [%s]',campanha,campanhacode),...
    'FontSize',12,'Fontweight','bold','interpreter','none');

h3=m_text(XL1,YL1,sprintf('%s [%s]',campanha,campanhacode), ...
    'Interpreter','none','HorizontalAlignment','center','FontSize',13,'FontWeight','bold');
h4=m_text(XL2,YL2,sprintf('%s [%s]',campanha,campanhacode), ...
    'Interpreter','none','HorizontalAlignment','center','FontSize',12);


CreaFigura(gcf,strcat('DatosCampanha',campanha),4)

%% Save file
clear BAT C h cbiso iiso
save('DatosCampanha')
