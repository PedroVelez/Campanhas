% ProcesaSBE911Btl - Post SEABIRD software procesing
%
%       Este programa visualiza  los perfiles de botellas del CTD sbe911
%       Despues calcula la densidad potencial referida 0dbar y lo pasa a formato .mat con
%       el mismo nombre que el fichero de entrada
%       Es el segundo paso despues de haber creado los .cnv con el software de SBE y
%       haber hecho el procesaSbe911.m para los .cnv
%
%		v1.0 15 Julio 2003 - Instituto Espanol de Oceanografia
%		v1.1 25 Junio 2004 - Instituto Espanol de Oceanografia
%		v1.2 25 Junio 2005 - Instituto Espanol de Oceanografia
%			Hacemos cambios para que no de errores cuando no hay ficheros de batimetria y/o costa
%
clear all;close all;clc
%-------------------------------------------------

load('../DatosCampanha')

presmin=8;
automatico=0;

%% Begin
r=1;
np=0;
Files=dir('*.btl');

while r==1;
    close all;
    np=np+1;
    keep r lat_min lat_max lon_min lon_max filebat filecosta filemalla presmin Files np automatico
    close all;
    file=Files(np).name;
    fprintf('>>>>>Procesando estacion %s  (%d/%d) \n',file,np,length(Files));
    [lati,long,nstat,gtime,depth,data,names]=Sbe911Btl_mat(file);
    numbtl=data.numbtl;
    prebtl=data.prebtl;
    salbtl=data.salbtl;
    if isfield(data,'oxybtl')
        oxybtl=data.oxybtl;
    else
        oxybtl=data.salbtl*NaN;
    end
    if isfield(data,'oxy2btl')
        oxy2btl=data.oxy2btl;
    else
        oxy2btl=data.salbtl*NaN;
    end
    if isfield(data,'turbtl')
        turbtl=data.turbtl;
    else
        turbtl=data.salbtl*NaN;
    end
    if isfield(data,'tur2btl')
        tur2btl=data.tur2btl;
    else
        tur2btl=data.salbtl*NaN;
    end
    if isfield(data,'flubtl')
        flubtl=data.flubtl;
    else
        flubtl=data.salbtl*NaN;
    end
    
    
    tembtl=data.tembtl;
    conbtl=data.conbtl;
    if max(numbtl)~=length(numbtl)
        numbtl=1:length(numbtl);
    end
    
    %----------
    %Figuras
    figure(1);%set(1,'position',[0.015 0.326667 0.5675 0.545])
    m_proj('Mercator','long',[lon_min lon_max],'lat',[lat_min lat_max]);hold on
    m_usercoast(filecosta,'patch',[.5 .5 .5]);
    m_grid('box','fancy','tickdir','in');
    title(sprintf('File: %s, estacion:%3.0d, \n Lat:%5.3f, Lon:%5.3f, depth(dbar): %4.1f \n Date:%s', ...
        file,nstat,lati,long,depth,datestr(datenum(gtime(1),gtime(2),gtime(3),gtime(4),gtime(5),gtime(6)))),'interpreter','none');
    m_line(long+360,lati,'marker','o','markersize',6,'color','b','MarkerFaceColor','b');
    m_text(long-0.15+360,lati+0.1,num2str(nstat),'Fontsize',8,'HorizontalAlignment','center','VerticalAlignment','top');
    
    
    figure(2);%set(gcf,'position',[0.59375 0.0566667 0.4025 0.82])
    %temperatura (perfil)
    subplot(1,3,1)
    plot(tembtl,prebtl,'bo-');hold on
    if exist('temp2btl','var')
        plot(temp2btl,prebtl,'r');
    end
    zoom on;grid on
    set(gca,'ydir','reverse')
    title('Temperatura')
    
    %salinidad (perfil)
    subplot(1,3,2)
    plot(salbtl,prebtl,'bo-');hold on
    if exist('salt2btl','var')
        plot(salt2btl,prebtl,'r');
    end
    zoom on;grid on
    set(gca,'ydir','reverse')
    title('Salinidad')
    
    %salinidad (perfil)
    subplot(1,3,3)
    plot(oxybtl,prebtl,'bo-');hold on
    zoom on;grid on
    set(gca,'ydir','reverse')
    title('Oxygen')
    
    
    save(strcat(file(1:length(file)-4),'btl'),'long','lati','nstat','gtime','depth','numbtl','prebtl','tembtl','salbtl','conbtl','oxybtl')
    if exist('sal2btl','var') && exist('tem2btl','var')
        save(strcat(file(1:length(file)-4),'btl'),'-append','tem2btl','sal2btl','con2btl')
    end
    
    if exist('turbtl','var')
        save(strcat(file(1:length(file)-4),'btl'),'-append','turbtl')
    end
    if exist('oxy2btl','var')
        save(strcat(file(1:length(file)-4),'btl'),'-append','oxy2btl')
    end
    if exist('tur2btl','var')
        save(strcat(file(1:length(file)-4),'btl'),'-append','tur2btl')
    end
    if exist('flubtl','var')
        save(strcat(file(1:length(file)-4),'btl'),'-append','flubtl')
    end
    
    
    fprintf('     >Datos almacenados en %s\n',strcat(file(1:length(file)-4),'btl','.mat'))
    
    if automatico==0 && np<length(Files)
        r=input('>>>>> Seguimos procesando ([1]/0)?');
        if r == 0
            break;
        end
        r=1;
    elseif size(Files,1)==np
        break;
    end
end