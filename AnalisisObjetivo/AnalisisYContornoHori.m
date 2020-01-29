% AnalisisYContorno - Script for objective data analysis and representation 
%		
%		Este programa hace analisis objetivo de una variable y lo representa
%		Tambien existe la posibilidad de calcular velocidades geostroficas y 
%		representarlas, partiendo del campo de altura dinámica calculado previamente.
%
% 1.0 14 Julio 2003 Instituto Español de Oceanografia (c)
% 1.2 26 Junio 2004 Instituto Español de Oceanografia (c)
%               Posibilidad de añadir velocidad geostrofica
% 1.3 26 Junio 2005 Instituto Español de Oceanografia (c)
%               Some additional checkings                        
%		Output:
%		hex,hexc- handels de excontour:
%		hmc	- handels de los puntos con los datos observados m_line: 
%
%		Input:
%			file=''; 		%Lectura de datos, variables: tems,sals,pres, dyns 
%			filemtopo='';		%Fichero donde se encuentra la topografia sobre la malla (mascara), hecho con Mtopo.m
%			filebat=''; 		%Fichero con la batimetria de la zona
%			filecosta='';		%Fichero con la costa de la zona
%			variable='Salinity';	%Variable to analyze
%			data=sals;		%Variable to analyze
%			profundidad=25;	    %Profundidad donde se escoge el analisis
%			automatico=0; 		%Variable para hacer contornos automaticos o no
%				ies=0.05; 	%Intevalo de contorno
%				escala=36.50:ies:38.30;		%Valores de contorno
%			etiquetas=1;  		%Variable para poner etiquetas en algunos contornos, (2) automatico
%				labelsetiquetas=[36.5:0.25:38.5]; %Etiquetas a poner
%			dtrend=0;   		%Variable para quitar la media al campo observado. Se usa para dyn
%
%			figuras=1; 		%Variable para hacer figuras: (0)Nada (1)jpeg (2)ps (3)emf (4)png
%			cl=jet;			%Escala de color
%			marcas=0;		%Variable para poner puntos con el valor observado
%			msize=5; 		%Tamaño de los puntos
%			bati=0; 		%Variable para añadir topografia
%				batiprof=[-1000 -100];	%Contornos de batimetria que se ponen
%			mascara=1; 		%Variable para usar las mascara de batimetria
%			titlef=1;		%Añade titulo con el nombre de la seccion
%			subtitlef=1;		%Pone un subtitulo con los campos maximos/minimos observados y analizados
%
%		Parametros para velocidades
%			fvector=1;		%Variable para poner vectores o no
%			munan=4;   		%Velocidades menores que munan no se representan
%			fac=0.006; 		%Factor de escala para los vectores
%			ns=3; 			%Number of velocityvectors to skip
%
%		Parametros del analisis
%			oberrormax=.3;		%Valor del analysis noise-to-signal ratio a partir del cual ya no considero el analisis
%			eta=0.05;		%Noise-to-signal ratio
%			corlen=0.2;%.3 		%Escala de Correlación
%			merror=1;		%Variable para poner o no el mapa de errores
%
%		Area geografica 
%			lat_min=37.75;
%			lat_max=40.5;
%			lon_min=0;
%			lon_max=5;
%
%		Area de analisis (Now read from DatosCampanha)
%			xnew=0.5:0.050:5.05;	%Mallado en x
%			ynew=37.8:0.050:40.40;	%Mallado en y
%
%----------------------------------------------------------------------
close all;clc

%----------------------------------------------------------------------
% Some checkings for the flags, if the flags are not defined the default 
% values are used
%----------------------------------------------------------------------
if exist('campanha','var')==0; campanha='';end
if exist('filebat','var')==0; filebat='';end
if exist('filecosta','var')==0; filecosta='';end
if exist('filemalla','var')==0; filemalla='';end
if exist('filemtopo','var')==0; filemtopo='';end

if exist('cl','var')==0; cl=jet; end
if exist('profundidad','var')==0;	profundidad=10;end
if exist('automatico','var')==0;	automatico=1;end
if exist('etiquetas','var')==0;	etiquetas=2;end
if exist('figuras','var')==0;	figuras=0;end
if exist('marcas','var')==0;	marcas=1;end
if exist('msize','var')==0;	msize=5;end
if exist('bati','var')==0;	bati=0;end
if exist('mascara','var')==0;	mascara=1;end
if exist('titlef','var')==0;	titlef=1;end
if exist('subtitlef','var')==0;	subtitlef=1;end
if exist('merror','var')==0;	merror=0;end
if exist('fvector','var')==0;   fvector=0;end
if exist('munan','var')==0;	    munan=3;end
if exist('dtrend','var')==0;    dtrend=0;end
if exist('fac','var')==0;	fac=0.006;end
if exist('ns','var')==0;	ns=1;end

if exist('oberrormax','var')==0  | exist('eta','var')==0 |exist('corlen','var')==0
    disp('>>>> One of the objective analysis parameteres has not been set')
    return
end    
if exist('xnew','var')==0  | exist('ynew','var')==0 
    disp('>>>> The analysis domain has not been set, xnew or ynew does not exist')
    return
end    

if automatico==0 & length(escala)==0
    automatico=1;
end
if automatico==0 & length(ies)==0
    automatico=1
end
if etiquetas==1 & length(labelsetiquetas)==0
    etiquetas=0;
end

%----------------------------------------------------------------------
%Begin
%----------------------------------------------------------------------
%Ask for the depth if necesary
if profundidad<0
    profundidad=input('>>>>> Depth? ');
end    
%Name of the section
if  length(findstr(file,'\'))==0
    nombresec=file;
else
    nombresec=file(max(findstr(file,'\'))+1:size(file,2));
end
disp(sprintf('>>>>> Analisis of %s at %03d dbar in %s',variable,profundidad,nombresec))

%Define the area of the analisis
[X,Y]=meshgrid(xnew,ynew);
%Lookfor the observations at the depth level
indice=find(pres == profundidad);
if length(indice)==0
    disp('Error: There is not data at the selected depth level')
    return
else
    zobs=data(indice,:)';
    yobs=lats';
    xobs=lons';
end    

%Quito aquellos lugares donde no hay datos
indice=find(isnan(zobs)-1); %posiciones donde hay datos
zobs=zobs(indice);
yobs=yobs(indice);
xobs=xobs(indice);
if length(zobs)==0
    disp('Error: There is not data')
    return
end    

%Quito la media, ya que el analisis es solo sobre residuos
m=mean(zobs);
zobs=zobs-m;

%Analisis propiamente dicho
disp(sprintf('     >Objective analysis with the following parameters\n     >eta:%5.2f, corlen:%5.2f',eta,corlen))
[psi,oberror] = objan(xobs,yobs,zobs,eta,corlen,xnew,ynew);
indice=find(oberror>oberrormax);
psi(indice)=NaN;

%Repongo la media
if exist('dtrend','var') ~= 0
    if dtrend==0 
        psi=psi+m;
        zobs=zobs+m;
    else
        disp('     >Quito la media a los campos medidos y analizados')
    end   
else
    psi=psi+m;
    zobs=zobs+m;
end
%---------------------
%Mapa de errores
%---------------------
if merror==1
    disp('     >Mapa de Errores')
    figure(1)
    m_proj('Mercator','long',[lon_min lon_max],'lat',[lat_min lat_max]);hold on
    [c,h]=m_contour(xnew,ynew,oberror,nice(extrem(oberror)));
    clabel(c,h); clear c h
    [lonsp,latsp]=m_ll2xy(lons,lats);
    plot(lonsp,latsp,'.');
    %Añadimos ficheros de costa
    %Añado la batimetría
    if bati == 1
        if length(filebat)>0
            load(filebat);
            m_extcontour(batylon,batylat,elevations,batiprof,'color',[0.75 0.75 0.75],'linewidth',1);
        end
    end
    if length(filecosta)>0
        m_usercoast(filecosta,'patch',[.5 .5 .5]);
    end      
    m_grid('box','fancy','tickdir','in');
end   
%---------------------
%Contorno del analisis
%---------------------
figure(2)
m_proj('Mercator','long',[lon_min lon_max],'lat',[lat_min lat_max]);hold on

%Elimino los datos donde no hay topografia. La matriz Mtopo se ha creado con la funcion Crea Mtopo.m
%teniendo en cuenta el grid que se usa, que debe de corresponder con el definido en DatosCampanha
if mascara == 1
    disp(sprintf('     >Elimino los datos donde no hay topografia'))
    if length(filemtopo)>0
        load(filemtopo);
        indice=find(H>(-profundidad-2));
        psi(indice)=NaN;
        clear indice H
    end
end
%La escala se pone automaticamente
if automatico == 1
    escala=nice(extrem([extrem(psi) extrem(zobs)]),4);
    disp('     >Escalado automatico')
end   
%Añado la batimetría
if bati == 1
    disp(sprintf('     >Poniedo batimetria en: %s',num2str(batiprof)))
    if length(filebat)>0
        load(filebat);
        m_extcontour(batylon,batylat,elevations,batiprof,'color',[0.75 0.75 0.75],'linewidth',1);
    end
end

%Contorno del analisis propiamente dicho
[hex,hexc]=m_extcontour(xnew,ynew,psi,escala,'fill');
caxis([min(escala) max(escala)])

%Pongo marcas donde hay estaciones con datos
if marcas == 1
    disp('     >Poninedo marcas donde hay datos')
    cl=jet;
    color=63*((zobs-min(escala))/(max(escala)-min(escala)))+1;
    color=63*((zobs-min(escala))/(max(escala)-min(escala)))+1;
    color(find(zobs>max(escala)))=64;
    color(find(zobs<min(escala)))=1;
    for i=1:size(xobs)
        hmc(i)=m_line(xobs(i),yobs(i),'color',cl(ceil(color(i)),:), ... 
            'marker','o','markersize',msize,'MarkerFaceColor',cl(ceil(color(i)),:), ...
            'MarkerEdgeColor','k','LineWidth',0.2);
    end
end

%Pongo etiquetas
if etiquetas == 1 | etiquetas == 2
    if etiquetas == 2 
        labelsetiquetas=nice(extrem([extrem(psi) extrem(zobs)]),1);
        labelsetiquetas=labelsetiquetas(1:2:end);
        disp(sprintf('     >Poniedo etiquetas automaticas en: %s',num2str(labelsetiquetas)))
    end
    [cs,h]=m_contour(xnew,ynew,psi,labelsetiquetas,'k');set(h,'LineWidth',1);
    ht=clabel(cs,h);set(ht,'HorizontalAlignment','center','VerticalAlignment','top','FontName','aerial','Rotation',0,'fontsize',6,'color','w');
end

%Añado Las velocidades geostroficas
if exist('fvector','var') == 1
    if fvector == 1
        if exist('dyns','var') ~= 0
            %AnalisisYContornoVec
            %Computing the geostrophic velocity from dynamic height
            %If the dyn height is in dyncm the geostrophic velocity is in cm/s
            disp(sprintf('     >Geostrophic velocity at %03d dbar in %s',profundidad,nombresec))
            %Escogo los valores en el nivel deseado
            indice=find(pres == profundidad);
            zobsdyn=10.*dyns(indice,:)'; %Factor 10 is to convert Vgeos to cm/s
            yobsdyn=lats';
            xobsdyn=lons';
            clear indice
            
            %Quito aquellos lugares donde no hay datos
            indice=find(isnan(zobsdyn)-1); %posiciones donde hay datos
            zobsdyn=zobsdyn(indice);
            yobsdyn=yobsdyn(indice);
            xobsdyn=xobsdyn(indice);
            clear indice
            
            %Quito la media, ya que el analisis es solo sobre residuos
            m=mean(zobsdyn);
            zobsdyn=zobsdyn-m;
            
            %Analisis de dyn propiamente dicho
            disp(sprintf('     >Realizando analisis de dyn con los siguientes parametros\n     >eta:%5.2f, corlen:%5.2f',eta,corlen))
            [psidyn,oberror] = objan(xobsdyn,yobsdyn,zobsdyn,eta,corlen,xnew,ynew);
            indice=find(oberror>oberrormax);
            psidyn(indice)=NaN;
            
            %Elimino los datos donde no hay topografia
            if mascara == 1
                disp(sprintf('     >Elimino los datos donde no hay topografia'))
                if length(filemtopo)>0
                    load(filemtopo);
                    indice=find(H>(-profundidad-10));
                    psidyn(indice)=NaN;
                    clear indice H
                end   
            end
            
            %Computing geostrophic velocity
            dx=sw_dist([ynew(1) ynew(1)],[xnew(1) xnew(2)],'km')*1000;
            dy=sw_dist([ynew(1) ynew(2)],[xnew(1) xnew(1)],'km')*1000;
            f=sw_f(mean(ynew)); %Coriolis factor
            [v,tu]=gradient(psidyn,dy,dx);
            u=-tu.*(10/f);
            v=v.*(10/f); %Si el campo de geopotencial estaba en cmdyns ahora obtengo cm/s
            mu=sqrt(u.^2+v.^2); %Modulo  del vector U
            Xm=X;Ym=Y;
            
            %Elimino unos datos para que el grafico sea mas claro
            %Velocities<munan=NaN
            %Field(oberror>oberrormax)=NaN
            indice=[find(oberror>oberrormax);find(mu<munan)]; 
            u(indice)=NaN; v(indice)=NaN;
            indice=[find((isnan(u)==1));find((isnan(v)==1))]; 
            Xm(indice)=NaN; Ym(indice)=NaN;
            %only one every ns values is represented
            um=u(1:ns:end,1:ns:end);vm=v(1:ns:end,1:ns:end);
            Xm=Xm(1:ns:end,1:ns:end);Ym=Ym(1:ns:end,1:ns:end);
            um=um(:);vm=vm(:);
            Xm=Xm(:);Ym=Ym(:);
            
            %Plotting the geostrophic velocity
            disp('     >Añadiendo velocidades geostroficas')
            hqv=m_quiver(Xm,Ym,um*fac,vm*fac,0,'k','filled');set(hqv,'linewidth',1);
            %Reference vector
            hqv=m_quiver(1.25,40.25,25.*fac,0.,0,'k','filled');set(hqv,'linewidth',1);
            ht=m_text(1.25,40.2,'25 cm/s','HorizontalAlignment','left','FontName','aerial','Rotation',0,'fontsize',6);
            
            disp(sprintf('     >Campo velocidad U analizado Min %5.2f Max %5.2f',min9(min9(u)),max9(max9(u))))
            disp(sprintf('     >Campo velocidad V analizado Min %5.2f Max %5.2f',min9(min9(v)),max9(max9(v))))   
        else   
            disp('>>>>> No existe la varible altura dinamica')
        end
    end
end

disp(sprintf('     >Campo observado Min %5.2f Max %5.2f',min(zobs),max(zobs)))
disp(sprintf('     >Campo analizado Min %5.2f Max %5.2f',min(min(psi)),max(max(psi))))

%Titles
if titlef == 1
    title(sprintf('%s at %s dbar %s',variable,num2str(profundidad),nombresec),'interpreter','none')
end   
if subtitlef == 1
    xlabel(sprintf('Observed field Min %5.2f Max %5.2f \nAnalyzed field Min %5.2f Max %5.2f',min(zobs),max(zobs),min(min(psi)),max(max(psi))),'HorizontalAlignment','center');
end   

%Colorbar
h=colorbar;set(h,'FontName','aerial','fontsize',8);
colormap(cl);caxis([min(escala) max(escala)])
%Incluyo la costa y el grid
if length(filecosta)>0
    m_usercoast(filecosta,'patch',[.5 .5 .5]);
end      
m_grid('box','fancy','tickdir','in','fontsize',10,'fontname','aerial');

%Making figures
creafigura(2,deblank(sprintf('%s_%s%03dm',nombresec,variable,profundidad)),figuras)
