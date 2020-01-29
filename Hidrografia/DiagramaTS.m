%Crea un DiagramaTS
clc;close all;clear

filedata='Ra1804';

temmax=23;temmin=2;
salmax=37;salmin=34.75;

Zona{1}.Est = [ 11: 24];    Zona{1}.tit='Norte';
Zona{2}.Est = [ 01: 10];    Zona{2}.tit='Canal Lanzarote';
Zona{3}.Est = [1101:1110];  Zona{3}.tit='La Graciosa';
Zona{4}.Est = [901:914];    Zona{4}.tit='Agadir';
Zona{5}.Est = [101:109];    Zona{5}.tit='Fuerteventura';
Zona{6}.Est = [502:507];    Zona{6}.tit='';

ZonaColor=['k' 'r' 'b' 'm' 'g' 'c'];
ZonaFlag =[1 1 1 1 1 1]; %Dibujo o no cada una de las zonas
Zonapsize=[2 2 2 2 2 2]; %Tamano de las marcas del TS

marcas=0; % Marcas sobre profundidaes seleccionadas
tipomarcap='o'; % Marcas sobre la TS
msize=Zonapsize*4;  % Tamano de las marcas
n=(500);        % Profundidades a la que se ponen las marcas

info=1;		    %Display por pantalla de datos del cast
choose =0;      %Pinto Estaciones concretas
figuras=[4 7]; 	%Variable para hacer figura (0)Nada (1)jpeg (2)ps (3)emf (4)png

x0t=36;
y0t=10;
dyt=1.1;

%% Inicio
load(filedata);
if exist('../DatosCampanha.mat','file') > 0
    DC=load('../DatosCampanha.mat');
else
    campanha='';
end

% En caso pintar estaciones concretas las pido
if choose,stchoose=input('Stations to plot?   ');end
ptmps = sw_ptmp(salts,temps,press,0);

%% Creo el Diagrama TS
figure(1)
% Creo la malla de densidad para el TS L?neas sg
tr = (temmin-1:0.5:temmax+1);		% set up meshgrid for
sr = (salmin-0.2:.1:salmax+0.2);	% plotting sig and ta
[S,T]=meshgrid(sr,tr);				% field in T/S plot
sg = sw_dens(S,T,zeros(size(T))) - 1000; 
[C,h]=contour(sr,tr,sg,[26.44 26.85 27.162 27.38 27.62 27.82 27.975],'-','LineWidth',.5,'color',[.65 .65 .65]);hold on
hc=clabel(C,h,[26.44 26.85 27.162 27.38 27.62 27.82],'fontsize',10);
titulo=DC.campanha;


%% Represento cada perfil
for jj=1:size(temps,2)
    pres=press;
    salt=salts(:,jj);
    temp=temps(:,jj);
    nstat=nstats(jj);
    if info
        fprintf('> Number of station %3.0d\n',nstats(jj))
        fprintf('  Depth of station: %4.0d\n',depths(jj));
        fprintf('  Date: %s\n',datestr(dates(jj)));
    end
    tpot=sw_ptmp(salt,temp/0.99976,pres,0); %El factor es para pasar de T90 a T68, que es la unidad que toman las SW
    if choose
        if any(nstat==stchoose)
            plot(salt,tpot,'k.');
        end
    else
        %Plot por zonas
        for iz=1:6
            if any(nstat==Zona{iz}.Est)
                if ZonaFlag(iz)==1,plot(salt,tpot,'o','MarkerEdgeColor',ZonaColor(iz),'MarkerFaceColor',ZonaColor(iz),'markersize',Zonapsize(iz));
                end
            end
        end %if
    end %if
    clear datevar gtime names par salt depth lati nsta pde sensors flu long oxy pres temp
end %for

if marcas
    %Represento cada perfil
    for jj=1:length(files)
        file=strcat(directorio,files(jj).name);
        load(file);
        if info
            fprintf('\n\nFile name %s\n',file)
            fprintf('Number of station %3.0d\n',nstat)
            fprintf('Coordinates (dec format): Lat: %6.4f   Lon: %6.4f\n',lati,long);
            fprintf('Depth of station: %4.0d\n',depth);
            fprintf('Date: %s\n\n',datestr(date));
        end %ifinfo
        if corregirsal
            if exist('sensors','var'),salt=repsal(salt,datevar,sensors);end %Aplicando correcci?n
        end %ifcorregirsal
        tpot=sw_ptmp(salt,temp/0.99976,pre,0); %El factor es para pasar de T90 a T68, que es la unidad que toman las SW
        if choose
            if any(noper==stchoose)
                plot(salt,tpot,'k.');
            end
        else
            for iz=1:6
                if any(noper==Zona{iz}.Est)
                    if ZonaFlag(iz)==1
                        if intersect(pre,n),plot(salt(intersect(pre,n)),tpot(intersect(pre,n)),tipomarcap,'color',ZonaColor{iz},'MarkerFaceColor',ZonaColor{iz},'marker',tipomarcap,'Markersize',msize(iz));end
                    end
                end %ifnoper
            end
        end %ifchoose
    end %for
end %ifmarcas
xlabel('Salinity','Fontsize',16)
ylabel('Potential Temperature','Fontsize',16)

if marcas==1
    title(sprintf('\\theta/S diagram %s Marks over %s db',titulo,num2str(n)),'Fontsize',14)
else
    title(sprintf('%s \\theta/S diagram %s',titulo),'Fontsize',14)
end
box on

%Titulo por zonas
if choose==0
    for iz=1:6
        if ZonaFlag(iz)==1
            text(x0t,y0t-(iz-1)*dyt      ,Zona{iz}.tit,'color',ZonaColor(iz),'FontSize',12,'FontWeight','bold','FontName','Verdana','BackgroundColor','w')
        end
    end
end

%Gamman
[C,h]=contour(salts,ptmps,gamas,[26.44 26.85 27.162 27.38 27.62 27.82 27.922 27.975],'k','linewidth',2,'color',[0.65 0.65 0.65]);
clabel(C,h,'fontsize',08,'LabelSpacing',600,'color',[0.65 0.65 0.65],'backgroundcolor','w')
grid on
CreaFigura(1,strcat(mfilename,'_',DC.campanha),figuras)
