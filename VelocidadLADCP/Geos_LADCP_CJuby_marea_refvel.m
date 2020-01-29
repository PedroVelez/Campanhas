% Lee velocidad gesotrofica de un par de estaciones, los valores LADCP de
% las dos estaciones, su media, grafica y compara el perficl de las dos
% velocidades
Limpia

DC=load('../DatosCampanha');
CodigoSeccion='CJuby';
sta=[108:-1:101]; %de oeste a este

angle=sta.*0-13.1;

CruiseDir=DC.dirdata(1:end-4);
vel_file=   strcat(CruiseDir,'/Analisis/Transporte/vel_',CodigoSeccion,'.mat');
marea_file= strcat(CruiseDir,'/Analisis/Marea/',DC.campanha,'_marea_',CodigoSeccion,'.mat');
QuitaMarea=1;
LADCP_file= strcat(CruiseDir,'/LADCP/Visbeck/profiles/');
output_file=strcat(CruiseDir,'/Analisis/VelocidadLADCP/refvel_',DC.campanha,'_',CodigoSeccion);

%%
%Numero de estaciones para LADCP. Quito las etsaciones del CTD que he
% eliminado
n_sta=length(sta);
load(vel_file)

% Marea
if QuitaMarea==1
    load(marea_file)
    % pasamos la velocidad a m/s
    TSU_L=TS_u*0.01;
    TSV_L=TS_v*0.01;
end

num_sta=size(stnvel,2)/2;

% En stnvel la variables impares son las profundidades y las pares las
% velocidades
[m,n]=size(stnvel);

depth=[];
vel_geos=[];
for ii=1:2:n-1
    depth=merge(depth,stnvel(:,ii));
    vel_geos=merge(vel_geos,stnvel(:,ii+1));
end

% LOs exactamente cero los pongo NaN. Es un problema del calculo
% geostrofico
index_depth_cero=find(depth==0);
depth(index_depth_cero)=NaN;
vel_geos(index_depth_cero)=NaN;

index_vel_cero=find(vel_geos==0);
depth(index_vel_cero)=NaN;
vel_geos(index_vel_cero)=NaN;

% LLamo al fichero LADCP
LADCP_vel_u=[];
LADCP_vel_v=[];
LADCP_z=[];
LADCP_vel_u_bottom=[];
LADCP_vel_v_bottom=[];
LADCP_z_bottom=[];
LADCP_lat=[];
LADCP_lon=[];

nombre=[];
for st =1:length(sta)
    if sta(st)<99
        flname=sprintf('%s%s_%03d.mat',LADCP_file,DC.campanhacode,sta(st));
    else
        flname=sprintf('%s%s_%03d.mat',LADCP_file,DC.campanhacode,sta(st));
    end
    fprintf('%s\n',flname)
    load(flname)
    nombre=strvcat(nombre,flname(end-6:end-4));
    LADCP_lon=merge(LADCP_lon,dr.lon);
    LADCP_lat=merge(LADCP_lat,dr.lat);
    LADCP_vel_u=merge(LADCP_vel_u,dr.u);
    LADCP_vel_v=merge(LADCP_vel_v,dr.v);
    LADCP_z=merge(LADCP_z,dr.z);
    if sta(st)==108 | sta(st)==102 | sta(st)==103 | sta(st)==101
        % Para las estaciones sin bot
        LADCP_vel_u_bottom=merge(LADCP_vel_u_bottom,NaN);
        LADCP_vel_v_bottom=merge(LADCP_vel_v_bottom,NaN);
        LADCP_z_bottom=merge(LADCP_z_bottom,NaN);
    else
        LADCP_vel_u_bottom=merge(LADCP_vel_u_bottom,dr.ubot);
        LADCP_vel_v_bottom=merge(LADCP_vel_v_bottom,dr.vbot);
        LADCP_z_bottom=merge(LADCP_z_bottom,dr.zbot);
    end
end

% Le quito la marea
if QuitaMarea==1
    [n,m]=size(LADCP_vel_u);
    TSU_1=repmat(TSU_L,n,1);
    TSV_1=repmat(TSV_L,n,1);
    LADCP_vel_u=LADCP_vel_u-TSU_1;
    LADCP_vel_v=LADCP_vel_v-TSV_1;
    % Y al bottom
    [n,m]=size(LADCP_vel_u_bottom);
    TSU_2=repmat(TSU_L,n,1);
    TSV_2=repmat(TSV_L,n,1);
    LADCP_vel_u_bottom=LADCP_vel_u_bottom-TSU_2;
    LADCP_vel_v_bottom=LADCP_vel_v_bottom-TSV_2;
end
% Paso a presion
LADCP_pres=sw_pres(LADCP_z,LADCP_lat);
LADCP_pres_bottom=sw_pres(LADCP_z_bottom,LADCP_lat);
refvel=NaN(1,length(sta)-1);

for ii=1:length(sta)-1
    LADCP_vel_u_mean(:,ii)=(LADCP_vel_u(:,ii)+LADCP_vel_u(:,ii+1))/2;
    LADCP_vel_v_mean(:,ii)=(LADCP_vel_v(:,ii)+LADCP_vel_v(:,ii+1))/2;
    LADCP_vel_u_bottom_mean(:,ii)=(LADCP_vel_u_bottom(:,ii)+LADCP_vel_u_bottom(:,ii+1))/2;
    LADCP_vel_v_bottom_mean(:,ii)=(LADCP_vel_v_bottom(:,ii)+LADCP_vel_v_bottom(:,ii+1))/2;
end

for ii=1:length(sta)-1
    LADCP_vel_u_cor(:,ii)=LADCP_vel_u_mean(:,ii).*cosd(angle(ii))+LADCP_vel_v_mean(:,ii).*sind(angle(ii));
    LADCP_vel_v_cor(:,ii)=-LADCP_vel_u_mean(:,ii).*sind(angle(ii))+LADCP_vel_v_mean(:,ii).*cosd(angle(ii));
    LADCP_vel_u_bottom_cor(:,ii)=LADCP_vel_u_bottom_mean(:,ii).*cosd(angle(ii))+LADCP_vel_v_bottom_mean(:,ii).*sind(angle(ii));
    LADCP_vel_v_bottom_cor(:,ii)=-LADCP_vel_u_bottom_mean(:,ii).*sind(angle(ii))+LADCP_vel_v_bottom_mean(:,ii).*cosd(angle(ii));
end

for ii=1:length(sta)-1
    LADCP_promedio=LADCP_vel_v_cor(:,ii);
    if     ii==1 ; mindepth=NaN; maxdepth=NaN;
    elseif ii==2 ; mindepth=NaN; maxdepth=NaN;
    elseif ii==3 ; mindepth=NaN; maxdepth=NaN;
    elseif ii==4 ; mindepth=NaN; maxdepth=NaN;
    elseif ii==5 ; mindepth=NaN; maxdepth=NaN;
    elseif ii==6 ; mindepth=NaN; maxdepth=NaN;
    elseif ii==7 ; mindepth=NaN; maxdepth=NaN;
    end
    
    ind_range=find(depth(:,ii)>=mindepth & depth(:,ii)<=maxdepth);
    selecc_vel=vel_geos(ind_range,ii);
    geovel_avg(ii)=nanmean(selecc_vel);
    
    ind_range2=find(LADCP_pres(:,ii)>=mindepth & LADCP_pres(:,ii)<=maxdepth);
    selecc_vel2=LADCP_promedio(ind_range2);
    u_profile(ii)=nanmean(selecc_vel2);
    
    refvel(ii)=u_profile(ii)-geovel_avg(ii);
    
    figure
    dd=plot(LADCP_promedio,-LADCP_z(:,ii),'b','linewidth',2);hold on;grid on
    plot(smooth(LADCP_promedio,100),-LADCP_z(:,ii),'b:','linewidth',1);hold on;grid on
    plot(LADCP_vel_v_bottom_cor(:,ii),-LADCP_z_bottom(:,ii),'r','linewidth',2)
    plot(LADCP_vel_v_bottom(:,ii+1),-LADCP_z_bottom(:,ii+1),'g','linewidth',2)
    plot(vel_geos(:,ii),-depth(:,ii),'-k','linewidth',2)
    plot(vel_geos(:,ii)+repmat(refvel(ii),length(depth(:,ii)),1),-depth(:,ii),'k','linewidth',4)
    xx=get(gca,'XLim');
    plot([xx(1) xx(2)],[-mindepth -mindepth],'linestyle','--','color',[.5 .5 .5],'linewidth',1.5)
    plot([xx(1) xx(2)],[-maxdepth -maxdepth],'linestyle','--','color',[.5 .5 .5],'linewidth',1.5)
    axis([-0.2 0.2 -inf inf])
    ylabel('Profundidad','fontsize',16)
    xlabel('Velocidad(m/s)','fontsize',16)
    
    title(sprintf('%s %03d-%03d; ii=%d',DC.campanhacode,sta(ii),sta(ii+1),ii),'fontsize',16);
    orient portrait
    CreaFigura(gcf,strcat(mfilename,sprintf('%03d-%03d',sta(ii),sta(ii+1))),7)
end

figure
plot(LADCP_lon(1:length(sta)-1),refvel,'*-b','markersize',12,'linewidth',1.5)
title('Refvel','fontsize',16)
ylim([-0.2 0.2])
grid on
ad1=add_scale('t',LADCP_lon(1:length(sta)),sta,'k');
set(ad1,'visible','off');

fprintf('Saving Results in %s\n',output_file)
save(output_file,'refvel','sta')
