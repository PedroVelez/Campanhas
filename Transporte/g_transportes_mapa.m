Limpia

Region='CanaryIslands';

Titulo='Raprocan1504';

DT=load('Raprocan1504_masa.mat');
DH=load('../../Datos/Ra1504.mat');

%Configuracion geografica
lat_min=26.5; lat_max=29.5;  lon_min=360-19;  lon_max=360-12;

load(strcat(GlobalSU.LibPath,'/Settings/DS',Region)); %Cargo Data Settings
load(GlobalDS.filebat)

%% Cre figuras
figure
m_proj('mercator','long',[lon_min lon_max],'lat',[lat_min lat_max]);hold on
m_contour(batylon,batylat,elevations,[-5000 -5000],'color',[0.25 0.25 0.25]);
m_contour(batylon,batylat,elevations,[-4000 -4000],'color',[0.35 0.35 0.35]);
m_contour(batylon,batylat,elevations,[-3000 -3000],'color',[0.45 0.45 0.45]);
m_contour(batylon,batylat,elevations,[-2000 -2000],'color',[0.55 0.55 0.45]);
m_contour(batylon,batylat,elevations,[-1000 -1000],'color',[0.65 0.65 0.65]);
m_contour(batylon,batylat,elevations,[ -500  -500],'color',[0.75 0.75 0.75]);
m_contour(batylon,batylat,elevations,[ -250  -250],'color',[0.85 0.85 0.85]);
m_usercoast(GlobalDS.filecoast,'patch',[.7 .6 .4,],'edgecolor',[.7 .6 .4,]);
m_grid('linestyle','none','fontsize',08,'fontname','verdana')

%Asigno variables
mass_trans=    [fliplr(DT.mass_trans_lanzarote)     ones(12,1)*0 fliplr(DT.mass_trans_norte)     -DT.mass_trans_oeste      -DT.mass_trans_sur];
mass_trans_ref=[fliplr(DT.mass_trans_lanzarote_ref) ones(12,1)*0 fliplr(DT.mass_trans_norte_ref) -DT.mass_trans_oeste_ref  -DT.mass_trans_sur_ref];

Lons=[DH.longs(1:end)];
Lats=[DH.latis(1:end)];


%% Calculo distancias y orintaciones de las estaciones
mod=1;
for ii=1:size(Lons,2)-1
    [di,phaseangle] = sw_dist([Lats(ii) Lats(ii+1)],[Lons(ii)+360 Lons(ii+1)+360]);
    %angle of line between stations with x axis (East).
    %Range of values are -180..+180. (E=0, N=90, S=-90)
    m_plot([Lons(ii)+360],[Lats(ii)],'ks')
    m_plot([Lons(ii)+360 Lons(ii+1)+360],[Lats(ii) Lats(ii+1)],'k-')
    dist(ii)=di;
    theta(ii)=phaseangle;
    LonsC(ii)=0.5*(Lons(ii)+Lons(ii+1));
    LatsC(ii)=0.5*(Lats(ii)+Lats(ii+1));
    m_plot(LonsC(ii)+360,LatsC(ii),'b.')
    mod=1;
end
m_plot(LonsC+360,LatsC,'bo')




%no Hago fliplr porque quiero hacer el acumulado desde el sur
%no pongo el signo negativo para hacer convenio geogr?fico pues al venir de
%calculos de modelo inverso, estaban en convenio caja.

mass_trans_sup_ref=sum(mass_trans_ref(1:4,:))/1e9;
mass_trans_int_ref=sum(mass_trans_ref(5:6,:))/1e9;
mass_trans_prof_ref=sum(mass_trans_ref(7:12,:))/1e9;
%Calculo el transporte acumulado
masscum_sup_ref=[cumsum(mass_trans_sup_ref')];
masscum_int_ref=[cumsum(mass_trans_int_ref')];
masscum_prof_ref=[cumsum(mass_trans_prof_ref')];
masscum_tot_ref=masscum_sup_ref+masscum_int_ref+masscum_prof_ref;

for ii=1:size(mass_trans_sup_ref,2)
    mod=mass_trans_sup_ref(ii);
    V_sup_ref(ii)=-mod.*cos(theta(ii)*pi/180);
    U_sup_ref(ii)=mod.*sin(theta(ii)*pi/180);
    %m_vec(4,LonsC(ii)+360,LatsC(ii),U(ii),V(ii));
    mod=masscum_sup_ref(ii);
    Vcum_sup_ref(ii)=-mod.*cos(theta(ii)*pi/180);
    Ucum_sup_ref(ii)=mod.*sin(theta(ii)*pi/180);
end%

m_vec(15,LonsC+360,LatsC,U_sup_ref,V_sup_ref, ...
	  'shaftwidth', 1, 'headlength',1 ,...
	  'EdgeColor','k');

%   m_vec(15,LonsC+360,LatsC,U_sup_ref,V_sup_ref, ...
% 	  'shaftwidth', 10, 'headlength', 0,...
% 	  'color',[0.75 0.75 0.75],'EdgeColor','k');
% fac=0.1;
% m_quiver(LonsC+360,LatsC,U_sup_ref*fac,V_sup_ref*fac,0);


% figure
% plot(LonsC(1:17),masscum_sup_ref(1:17),'k'); hold on
% bar(LonsC(1:17),mass_trans_sup_ref(1:17),1)
CreaFigura(1,strcat(mfilename,campanha),7)