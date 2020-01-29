% CreaSecCtd
%
%		Programa para agrupar los casts de CTD que se deseen en formato matlab.
%		Las secciones asi creadas luego seran usadas para los calculos posteriores
%
%
%		v1.0 25 Junio 2003 - Instituto Espa?ol de Oceanografia
%		v1.1 25 Junio 2004 - Instituto Espa?ol de Oceanografia
%		v1.2 25 Junio 2005 - Instituto Espa?ol de Oceanografia
%			Hacemos cambios para que no de errores cuando no hay ficheros de batimetria y/o costa
close all;clear all;clc
%-------------------------------------------------

load('../DatosCampanha')
if exist('campanha','var')==0; campanha='';end
if exist('filebat','var')==0; filebat='';end
if exist('filecosta','var')==0; filecosta='';end
if exist('filemalla','var')==0; filemalla='';end
if exist('filemtopo','var')==0; filemtopo='';end

%-------------------------------------------------
%% Begin
%-------------------------------------------------
stns=input('>>>> Stations for the sections? ([1 2 ...])?     ') ;
dates=[];depths=[];nstats=[];latis=[];longs=[];
numsbtl=[];presbtl=[];consbtl=[];temsbtl=[];salsbtl=[];oxysbtl=[];oxy2sbtl=[];tursbtl=[];tur2sbtl=[];flusbtl=[];

i1=0;
for jj=1:length(stns)
    fileBtl=strcat('./Mat/',campanhacode,sprintf('_%02dbtl',stns(jj)),'.mat');
    if exist(fileBtl,'file') >0
        i1=i1+1;
        fprintf('     > Reading %s \n',num2str(stns(jj)))
        load(fileBtl);
        en=0;
        latis(i1)=lati;
        longs(i1)=long;
        depths(i1)=depth;
        nstats(i1)=nstat;
        dates(i1)=datenum(gtime(1),gtime(2),gtime(3),gtime(4),gtime(5),gtime(6));
        numsbtl=CatWithNans(numsbtl,i1,numbtl(:));
        
        % Concatenamos los valores del lanzado en sus posiciones (meollo)
        presbtl=CatWithNans(presbtl,i1,prebtl(:));
        temsbtl=CatWithNans(temsbtl,i1,tembtl(:));
        salsbtl=CatWithNans(salsbtl,i1,salbtl(:));
        consbtl=CatWithNans(consbtl,i1,conbtl(:));
        oxysbtl=CatWithNans(oxysbtl,i1,oxybtl(:));
        oxy2sbtl=CatWithNans(oxy2sbtl,i1,oxy2btl(:));
        tursbtl=CatWithNans(tursbtl,i1,turbtl(:));
        tur2sbtl=CatWithNans(tur2sbtl,i1,tur2btl(:));
        flusbtl=CatWithNans(flusbtl,i1,flubtl(:));
    end
end

figure(1)
m_proj('Mercator','long',[lon_min lon_max],'lat',[lat_min lat_max]);
if ~isempty(filebat)
    load(filebat)
    m_contour(batylon,batylat,elevations,[-6000 -6000],'color',[0.15 0.15 0.15]);hold on
    m_contour(batylon,batylat,elevations,[-5000 -5000],'color',[0.25 0.25 0.25]);
    m_contour(batylon,batylat,elevations,[-4000 -4000],'color',[0.35 0.35 0.35]);
    m_contour(batylon,batylat,elevations,[-3000 -3000],'color',[0.45 0.45 0.45]);
    m_contour(batylon,batylat,elevations,[-2000 -2000],'color',[0.55 0.55 0.45]);
    m_contour(batylon,batylat,elevations,[-1000 -1000],'color',[0.65 0.65 0.65]);
    m_contour(batylon,batylat,elevations,[ -500  -500],'color',[0.75 0.75 0.75]);
    m_contour(batylon,batylat,elevations,[ -250  -250],'color',[0.85 0.85 0.85]);
end
if ~isempty(filecosta)
    m_usercoast(filecosta,'patch',[.7 .6 .4,],'edgecolor',[.7 .6 .4,]);
end
for ins=1:length(nstats)
    m_line(longs(ins),latis(ins),'marker','o','markersize',2,'color','b','MarkerFaceColor','b');
    m_text(longs(ins)+360,latis(ins),num2str(nstats(ins)),'Fontsize',9,'HorizontalAlignment','center','VerticalAlignment','top');
end %for
m_grid('linestyle','none','fontsize',08,'fontname','tahona')



nombre=input('>>>>> Name of the section?','s');
title(strcat(campanha,' seccion:   ',nombre,'.mat'))
save(strcat(char(nombre),'btl'),'nstats','latis','longs','dates','depths','numbtl','prebtl','temsbtl','salsbtl','consbtl')

print(1,'-djpeg',strcat('Mapa_',nombre))
respuesta=0;
respuesta=input('>>>>> Do you want to save the data in an excel file? (1/[0])?');
if respuesta == 1
    %     xlswrite(strcat(char(nombre),'btl'),nstats, 'nstats')
    %     xlswrite(strcat(char(nombre),'btl'),numsbtl,'numsbtl')
    %     xlswrite(strcat(char(nombre),'btl'),presbtl,'presbtl')
    %     xlswrite(strcat(char(nombre),'btl'),saltsbtl,'saltsbtl')
    %     xlswrite(strcat(char(nombre),'btl'),tempsbtl,'tempsbtl')
    %     xlswrite(strcat(char(nombre),'btl'),tempsbtl,'condsbtl')
    %     if exist('salt2btl','var')
    %         xlswrite(strcat(char(nombre),'btl'),saltsbtl,'salt2sbtl')
    %         xlswrite(strcat(char(nombre),'btl'),tempsbtl,'temp2sbtl')
    %         xlswrite(strcat(char(nombre),'btl'),tempsbtl,'cond2sbtl')
    %     end
    t1=(nstats'*ones(24,1)')';
    numsbtl2=[1:1:24]'*ones(size(nstats));
    fid=fopen(strcat(strcat(char(nombre),'btl'),'.csv'),'w');
    fprintf(fid,'%%Station;   Bottle;    Latitude;    Longitude;         Date;          Time;          PrDM;        Sal00;        T090C;         C0S/m;     Sbeox0ML/L;     Sbox0Mm/Kg;     CTurbWETntu0;     SeaTurbMtr;      FlECO-AFL \n');
    
    
    for i2=1:size(nstats,2)
        for i1=1:size(numsbtl,1)
            fprintf(fid,'%03d;         %02d;      %7.2f;      %7.2f;      %s;      %s;      %7.2f;      %7.4f;      %7.4f;      %7.4f;      %7.4f;      %7.4f;      %7.4f;      %7.4f;      %7.4f \n',nstats(i2),numsbtl(i1,i2),latis(i2),longs(i2),datestr(dates(i2),23),datestr(dates(i2),13), presbtl(i1,i2),salsbtl(i1,i2),temsbtl(i1,i2),consbtl(i1,i2),oxysbtl(i1,i2),oxy2sbtl(i1,i2),tursbtl(i1,i2),tur2sbtl(i1,i2),flusbtl(i1,i2));
        end
    end
    fclose(fid);
end