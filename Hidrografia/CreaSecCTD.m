% CreaSecCtd
%		Programa para agrupar los casts de CTD que se deseen en formato matlab.
%		Las secciones asi creadas luego seran usadas para los calculos posteriores
%		v1.0 25 Junio 2003 - Instituto Espanol de Oceanografia
%		v1.1 25 Junio 2004 - Instituto Espanol de Oceanografia
%		v1.2 25 Junio 2005 - Instituto Espanol de Oceanografia
%			Hacemos cambios para que no de errores cuando no hay ficheros de batimetria y/o costa
close all;clear all;clc
%-------------------------------------------------


%% Begin
if exist('../DatosCampanha.mat','file') > 0
    DC=load('../DatosCampanha.mat');
    campanha=DC.campanha;
else
    campanha='';
    DC.campanhacode='';
end

stns=input('>>>> Estaciones para crear la secci?n ([1 2 ...])?     ') ;

dates=[];latis=[];longs=[];press=[];temps=[];salts=[];temp2s=[];salt2s=[];sgths=[];
flus=[];depths=[];nstats=[];ncasts=[];gamas=[];oxyis=[];

fprintf('     > Procesing %d stations \n',length(stns))
fprintf('     > stations ')

for ins=1:length(stns)
    fprintf('%4d ',stns(ins))
    file=strcat('../../CTD/Mat/',DC.campanhacode,'_',sprintf('%03i',stns(ins)));
    D=load(file);
    if D.PRES(1)~=2
        depth0=min(D.PRES);
        D.PRES=(2:2:max(D.PRES))';
        temp=cat(1,(NaN.*ones(1,(depth0-2)/2))',D.TEMP);
        salt=cat(1,(NaN.*ones(1,(depth0-2)/2))',D.SALT);
        if isfield(D,'oxyi')
            oxyi=cat(1,(NaN.*ones(1,(depth0-2)/2))',D.oxyi);
        end
        if isfield(D,'TEMP2')
            temp2=cat(1,(NaN.*ones(1,(depth0-2)/2))',D.TEMP2);
        end
        if isfield(D,'SALT2')
            salt2=cat(1,(NaN.*ones(1,(depth0-2)/2))',D.SALT2);
        end
    end
    if isempty(press) || (max(D.PRES)>=max(press))
        press=(2:2:max(D.PRES))';
    end
    latis(ins)=D.LATI;
    longs(ins)=D.LONG-360;
    depths(ins)=D.depth;
    nstats(ins)=D.nstat;
    namess{ins}=D.names;
    dates(ins)=datenum(D.gtime);
    % Concatenamos los valores del lanzado en sus posiciones (meollo)
    temps=CatWithNans(temps,ins,temp);
    salts=CatWithNans(salts,ins,salt);
    if exist('temp2','var')
        temp2s=CatWithNans(temp2s,ins,temp2);
    end
    if exist('temp2','var')
        salt2s=CatWithNans(salt2s,ins,salt2);
    end
    if exist('oxyi','var')
        oxyis=CatWithNans(oxyis,ins,oxyi);
    end
end

clear temp salt oxyi temp2 salt2
%% interpola para eliminar valores ausentes intermedios en los perfiles
np=size(temps,2);
presi=press;
%Creo la matriz de presiones, igual en dimensiones a la de temps y salts
press=press*ones(1,size(temps,2));
press(isnan(temps))=NaN;
press(isnan(salts))=NaN;
ndpress=temps*NaN;
ndtemps=temps*NaN;
ndsalts=temps*NaN;
ndtemp2s=temps*NaN;
ndsalt2s=temps*NaN;

fprintf('\n     > Interpolating ...')
for ip=1:np
    fprintf(' %02d ',ip)
    
    
    mpres=max(press(:,ip));
    
    temp=temps(:,ip);
    pres=press(:,ip);
    in=find(isnan(temp)==0 & isnan(pres)==0);
    temp=temp(in); pres=pres(in);
    ntemp=interp1(pres,temp,presi,'linear',NaN)';
    ndtemps(:,ip)=ntemp;
    
    salt=salts(:,ip);
    pres=press(:,ip);
    in=find(isnan(salt)==0 & isnan(pres)==0);
    salt=salt(in); pres=pres(in);
    nsalt=interp1(pres,salt,presi,'linear',NaN)';
    ndsalts(:,ip)=nsalt;
    

    temp2=temp2s(:,ip);
    pres=press(:,ip);
    in=find(isnan(temp2)==0 & isnan(pres)==0);
    temp2=temp2(in); pres=pres(in);
    ntemp2=interp1(pres,temp2,presi,'linear',NaN)';
    ndtemp2s(:,ip)=ntemp2;
    
    
    salt2=salt2s(:,ip);
    pres=press(:,ip);
    in=find(isnan(salt2)==0 & isnan(pres)==0);
    salt2=salt2(in); pres=pres(in);
    nsalt2=interp1(pres,salt2,presi,'linear',NaN)';
    ndsalt2s(:,ip)=nsalt2;
    
    if exist('oxyi','var')
        oxyi=oxyis(:,ip);
        pres=press(:,ip);
        in=find(isnan(oxyi)==0 & isnan(pres)==0);
        oxyi=oxyi(in); pres=pres(in);
        noxyi=interp1(pres,oxyi,presi,'linear',NaN)';
        ndoxyis(:,ip)=noxyi;
    end
end
fprintf(' .\n ',ip)
salts=ndsalts;
temps=ndtemps;
press=ndpress;
if exist('oxyi','var')
    oxyis=ndoxyis;
end
salt2s=ndsalt2s;
temp2s=ndtemp2s;
sgths=sw_pden(salts,temps,press,0)-1000;
press=presi;
gamas = gamma_n(salts,temps/0.99976,press*ones(1,size(temps,2)),longs,latis);
clear ndsalt2s ndtemp2s ndsalts ndtemps ndpress pressold

%-------------------------------------------------
figure(1)
m_proj('Mercator','long',[DC.lon_min DC.lon_max],'lat',[DC.lat_min DC.lat_max]);
if ~isempty(DC.filebat)
    load(DC.filebat)
    m_contour(batylon,batylat,elevations,[-6000 -6000],'color',[0.15 0.15 0.15]);hold on
    m_contour(batylon,batylat,elevations,[-5000 -5000],'color',[0.25 0.25 0.25]);
    m_contour(batylon,batylat,elevations,[-4000 -4000],'color',[0.35 0.35 0.35]);
    m_contour(batylon,batylat,elevations,[-3000 -3000],'color',[0.45 0.45 0.45]);
    m_contour(batylon,batylat,elevations,[-2000 -2000],'color',[0.55 0.55 0.45]);
    m_contour(batylon,batylat,elevations,[-1000 -1000],'color',[0.65 0.65 0.65]);
    m_contour(batylon,batylat,elevations,[ -500  -500],'color',[0.75 0.75 0.75]);
    m_contour(batylon,batylat,elevations,[ -250  -250],'color',[0.85 0.85 0.85]);
end
if ~isempty(DC.filecosta)
    m_usercoast(DC.filecosta,'patch',[.7 .6 .4,],'edgecolor',[.7 .6 .4,]);
end
for ins=1:length(nstats)
    m_line(longs(ins)+360,latis(ins),'marker','o','markersize',2,'color','b','MarkerFaceColor','b');
    m_text(longs(ins)+360,latis(ins),num2str(nstats(ins)),'Fontsize',9,'HorizontalAlignment','center','VerticalAlignment','top');
end %for
m_grid('linestyle','none','fontsize',08,'fontname','tahona')
nombre=input('>>>>> Nombre de la seccion:','s');
title(strcat(DC.campanha,', seccion:   ',nombre,'.mat'),'interpreter','none')

save(char(nombre),'campanha','nstats','depths','latis','longs','dates','press','temps','salts','sgths','gamas')
if exist('temp2s','var')
    save(char(nombre),'-append','temp2s','salt2s')
end
if exist('oxyis','var')
    save(char(nombre),'-append','oxyis')
end

print(1,'-djpeg',strcat(nombre,'_Mapa'))

r=input('>>>>> Paso los datos a un fichero excel (1/[0])?');
if r == 1
    fid=fopen(strcat(strcat(char(nombre)),'.csv'),'w');
    fprintf(fid,'%Station; Latitude;    Longitude         Date;          Time;          PrDM;        Sal00;        T090C;  \n');
    for iest=1:size(nstats,2)
        for iniv=1:size(temps,1)
            fprintf(fid,'%03d;     %7.2f;      %7.2f;      %s;      %s;      %7.2f;      %7.4f;      %7.4f;\n',nstats(iest),latis(iest),longs(iest),datestr(dates(iest),23),datestr(dates(iest),13), press(iniv),salts(iniv,iest),temps(iniv,iest));
        end
    end
    fclose(fid);
end

