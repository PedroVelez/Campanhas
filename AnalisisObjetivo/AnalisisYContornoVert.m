% SeccionVar - Script for vertical section representation 
%		
%		Este programa representa una seccion vertical, partiendo de una
%		seccion creada
%
% Input:
%   file='';load(file)  - Lectura de datos. Si file='' pregunta por el nombre del fichero
%   variable='';        - Variables: tems,sals,pres, dyns
%   data=sals;
%	tiposeccion=1; 	    - 1 seccion WE (varx=lon),2 seccion NS (varx=lat), 3 Classical         
%	automatico=0; 		- Variable para hacer contornos automaticos o no
%		ies=0.01; 		- Intevalo de contorno
%		escala=38.40:ies:38.6; %Valores de contorno
%	etiquetas=1;  		- Variable para poner etiquetas en algunos contornos, 2 para automatico
%	labelsetiquetas=[38.2,38.45,38.50,38.55,38.60]; %Etiquetas a poner
%	lineasad=0;  		- Variable para lineas adicionales
%	    labelslineasad=[38.25,38.35,38.45,38.55];
%	ejespres=[-200 -650];- Limite de los ejes en si tiene tres valores, los ejes se dividen en dos
%                          [0 -100 -350]
%	figuras=4;    		- Variable para hacer figura en 1 .jpeg, 2 .emf, 3 .ps, 4 .png
%	bati=0; 			- Variable para añadir topografia a la seccion
%	mposicion=1;        - Mapa con la posicion de las estaciones
%	cl=jet;       		- Colores para la barra
%
% 1.0 14 Julio 2003 Instituto Español de Oceanografia (c) 
% 1.2 26 Junio 2004 Instituto Español de Oceanografia (c)
% 1.3 25 Junio 2005 Instituto Español de Oceanografia (c)
%									
%----------------------------------------------------------------------

%----------------------------------------------------------------------
% Some checkings for the flags, if the flags are not defined the default 
% values are used
%----------------------------------------------------------------------
if exist('campanha','var')==0; campanha='';end
if exist('filebat','var')==0; filebat='';end
if exist('automatico','var')==0;	automatico=1;end
if exist('etiquetas','var')==0;	etiquetas=0;end
if exist('lineasad','var')==0;	lineasad=0;end
if exist('figuras','var')==0;	figuras=0;end
if exist('bati','var')==0;	bati=1;end
if exist('mposicion','var')==0;	mposicion=1;end
if exist('cl','var')==0;	cl=jet;end
if exist('tiposeccion','var')==0; tiposeccion=1;end
if tiposeccion<0
    tiposeccion=input('>>>>> Tiposeccion(1/2/3)?');
    if exist('tiposeccion','var')==0 | tiposeccion>3
        disp('>>>> The section type is not defined')   
        return
    end    
end    
if exist('ejespres','var')==0; ejespres=[0 -max(pres(:))];end    
if length(escala)==0;automatico=1;end
if length(ies)==0;automatico=1;end
if length(labelsetiquetas)==0;etiquetas=0;end
if length(labelslineasad)==0;lineasad=0;end
%----------------------------------------------------------------------
% Begin
%----------------------------------------------------------------------
% Creo la matriz de presiones, latitutdes, longitudes 
lons=ones(size(data,1),1)*lons;
pres=pres*ones(1,size(data,2));
lats=ones(size(data,1),1)*lats;
stats=ones(size(data,1),1)*[1:1:length(nopers)];

%Nombre de la seccion
if length(findstr(file,'\'))>0
    nombresec=file(max(findstr(file,'\'))+1:size(file,2));
else
    nombresec=file;
end

%Ejes segun el tipo de seccion
if tiposeccion==1 
    varx=lons;
    disp('>>>>> WE Longitudinal section')
elseif tiposeccion==2
    varx=lats;
    disp('>>>>> NS Meridional section')
elseif tiposeccion==3
    varx=stats;
    bati=0;
    disp('>>>>> Classical section ')
end   
%Ordeno la seccion para que varx sea creciente
[varx2,indicen]=sort(varx,2);
n=size(data,1);
for j = 1:n, 
    varx(j,:)=varx(j,indicen(j,:));
    data(j,:)=data(j,indicen(j,:));
    lons(j,:)=lons(j,indicen(j,:));
    lats(j,:)=lats(j,indicen(j,:));
end
nopers=nopers(indicen(1,:));

%Creo batimetria
if bati == 1 
    if length(filebat)>0
        load(filebat)
        lonsf=[];
        latsf=[];
        %Creo la batimetria
        if tiposeccion==1 
            for i=1:size(lons,2)-1
                lonsf= [lonsf lons(1,i):(lons(1,i+1)-lons(1,i))/20:lons(1,i+1)];
            end   
            latsf=lonsf*0+mean(lats(1,:));        
        elseif tiposeccion==2
            for i=1:size(lons,2)-1
                latsf= [latsf lats(1,i):(lats(1,i+1)-lats(1,i))/20:lats(1,i+1)];
            end   
            lonsf=latsf*0+mean(lons(1,:));   
        end   
        for i=1:size(latsf,2)
            batprof(i)=elevations(locate(batylat,latsf(i)),locate(batylon,lonsf(i)));
        end   
    end    
end    

%Creo etiquetas
if etiquetas == 2 
    labelsetiquetas=nice(extrem(data),1);
    labelsetiquetas=labelsetiquetas(1:2:end);
    disp(sprintf('     >Poniedo etiquetas automaticas en: %s',num2str(labelsetiquetas)))
end

%---------------
%Mapa de poscion
%---------------
if mposicion==1
    figure(1)
    m_proj('Mercator','long',[lon_min lon_max],'lat',[lat_min lat_max]);
    m_grid('box','fancy','tickdir','in');
    if length(filebat)>0
        load(filebat)
        m_extcontour(batylon,batylat,elevations,[-2000 -1000 -200 -100],'k','linewidth',1);hold on
    end    
    if length(filecosta)>0
        m_usercoast(filecosta,'patch',[.5 .5 .5]);
    end    
    for jj=1:length(nopers)
        m_line(lons(1,jj),lats(1,jj),'marker','o','markersize',3,'color','b','MarkerFaceColor','b');
        m_text(lons(1,jj),lats(1,jj)+0.15,num2str(nopers(jj)),'Fontsize',8,'HorizontalAlignment','center','VerticalAlignment','top');
    end %for
    m_grid('box','fancy','tickdir','in');
    title(sprintf('%s, %s',variable,nombresec),'Interpreter','none')
end

%---------------
%Vertical Section
%---------------
figure(2)
disp(sprintf('     >%s Max %5.2f Min %5.2f',variable,nanmax(nanmax(data)),nanmin(nanmin(data))))

%La escala se pone automaticamente
if automatico == 1
    escala=nice(extrem(data),4);
    disp('     >Escalado automatico')
end   
if length(ejespres)==2
    extcontour(varx,-pres,data,'fill',escala);hold on;
    axis([-inf inf ejespres(2) ejespres(1)])
    %Añado la batimetrí
    if bati == 1 
        if length(filebat)>0
            if tiposeccion==1 
                plot(lonsf,batprof);
            elseif tiposeccion==2
                plot(latsf,batprof)
            end   
        end    
    end
    
    %Añado una marca donde estan las estaciones
    plot(varx(1,:),ones(1,size(varx,2))*(min(ejespres)),'bs','MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','b');
    for i=1:size(nopers,2)
        text(varx(1,i),(min(ejespres)+15),num2str(nopers(i)),'VerticalAlignment','top','HorizontalAlignment','center','FontSize',8);
    end
    
    %Pongo etiquetas
    if etiquetas == 1 | etiquetas == 2
        if length(labelsetiquetas)>0
            [cs,h]=contour(varx,-pres,data,labelsetiquetas,'k');set(h,'LineWidth',1);
            ht=clabel(cs,h,'labelspacing',40);set(ht,'HorizontalAlignment','center','VerticalAlignment','top','FontName','aerial','Rotation',0,'fontsize',6,'rotation',0);
        end    
    end
    
    %Pongo lineas adicionales
    if lineasad == 1
        if length(labelslineasad)>0
            [cs,h]=contour(varx,-pres,data,labelslineasad,'k');set(h,'LineWidth',1.25);
        end
    end
    ylabel('Pressure [db]')
    title(sprintf('%s at %s',variable,nombresec),'Interpreter','none')
    h=colorbar;set(h,'FontName','aerial','fontsize',8);
    colormap(cl);caxis([min(escala) max(escala)])
    
elseif length(ejespres)==3
    %----------
    %Superficial section
    %----------
    ax1pos=[0.10 0.70 0.77 0.22];    
    subplot('position',ax1pos);ax1=gca;
    extcontour(varx,-pres,data,escala,'filled');hold on;set(ax1,'Xtick',[]);
    axis([-inf inf ejespres(2) ejespres(1)])
    %Etiquetas
    if etiquetas == 1 | etiquetas == 2
        if length(labelsetiquetas)>0
            [cs,h]=contour(varx,-pres,data,labelsetiquetas,'k');set(h,'LineWidth',1);
            ht=clabel(cs,h,'labelspacing',40);set(ht,'HorizontalAlignment','center','VerticalAlignment','top','FontName','aerial','Rotation',0,'fontsize',6,'rotation',0);
        end    
    end
    %Pongo lineas adicionales
    if lineasad == 1
        if length(labelslineasad)>0
            [cs,h]=contour(varx,-pres,data,labelslineasad,'k');set(h,'LineWidth',1.25);
        end
    end
    %Añado una marca donde estan las estaciones
    plot(varx(1,:),ones(1,size(varx,2))*ejespres(2),'bs','MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','b');
    for i=1:size(nopers,2)
        text(varx(1,i),1.08*ejespres(2),num2str(nopers(i)),'VerticalAlignment','top','HorizontalAlignment','center','FontSize',6);
    end
    caxis([min(escala) max(escala)])
    ylabel('Pressure [db]')
    title(sprintf('%s at %s',variable,nombresec),'Interpreter','none')
    
    %--------
    %deep section
    %--------
    ax2pos=[0.10 0.11 0.77 0.55];
    subplot('position',ax2pos);ax2=gca;

    extcontour(varx,-pres,data,escala,'filled');hold on;set(ax1,'Xtick',[]);
    axis([-inf inf ejespres(3) ejespres(2)])
    %Añado una marca donde estan las estaciones
    plot(varx(1,:),ones(1,size(varx,2))*ejespres(2),'bs','MarkerSize',3,'MarkerEdgeColor','k','MarkerFaceColor','b');
    %Etiquetas
    if etiquetas == 1 | etiquetas == 2
        if length(labelsetiquetas)>0
            [cs,h]=contour(varx,-pres,data,labelsetiquetas,'k');set(h,'LineWidth',1);
            ht=clabel(cs,h,'labelspacing',40);set(ht,'HorizontalAlignment','center','VerticalAlignment','top','FontName','aerial','Rotation',0,'fontsize',6,'rotation',0);
        end    
    end
    %Pongo lineas adicionales
    if lineasad == 1
        if length(labelslineasad)>0
            [cs,h]=contour(varx,-pres,data,labelslineasad,'k');set(h,'LineWidth',1.25);
        end
    end
    %Añado la batimetrí
    if bati == 1 
        if length(filebat)>0
            if tiposeccion==1 
                plot(lonsf,batprof);
            elseif tiposeccion==2
                plot(latsf,batprof)
            end   
        end    
    end
    
    ax3=colorbar;set(ax3,'FontName','aerial','fontsize',8);
    colormap(cl);    
    caxis([min(escala) max(escala)])
    set(ax3,'position',[0.89,ax2pos(2),0.06,(ax1pos(4)+ax1pos(2))-ax2pos(2)])
    set(ax2,'position',ax2pos)
    ylabel('Pressure [db]')
end    
if tiposeccion==1 
    xlabel(sprintf('Longitude \nMax %5.2f Min %5.2f',nanmax(nanmax(data)),nanmin(nanmin(data))),'Interpreter','none')
elseif tiposeccion==2
    xlabel(sprintf('Latitude \nMax %5.2f Min %5.2f',nanmax(nanmax(data)),nanmin(nanmin(data))),'Interpreter','none')
end   

orient portrait
creafigura(2,sprintf('%s_%s%s',variable,nombresec),figuras)