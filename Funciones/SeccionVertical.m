function [ax1,ax2,ax3]=SeccionVertical(Data,Opciones);
%Data.xvar
%Data.yvar
%Data.zvar
%Opciones.tipo='Z';   %Tipo de seccion: 'Z' Zonal [Defecto], 'M' meridional
%                      'T' temporal
%Opciones.xrange=0;   %Rango ejeX - 0 para autom?tico o vector [-21 -12.5];
%Opciones.yrange1=0;  %Rango ejeY Superior- 0 para autom?tico o vector [-400 0];
%Opciones.yrange2=0;  %Rango ejeY inferior- 0 para autom?tico o vector [-1500 -400];
%Opciones.zrange=0;   %Rango de Z  - 0 para automatico o vector
%Opciones.ziso=0;     %isolineas a dibujary etiquetar - 0 para automatico o vector
%Opciones.ziso2=0;    %isolineas a etiquetar en --
%Opciones.stations=   %No pinta estaciones
%Opciones.batymetry=  % elevations values (negative
%                     % filename with variables batylon, elevations
%Opciones.colormap=   %

%% Comprobaciones
fprintf('     >>>>> SeccionVertical\n');
%% Data
xvar=Data.xvar;
yvar=Data.yvar;
zvar=Data.zvar;
if size(xvar,1)==1
    xvar=xvar';
end
%% Opciones

if strcmp(Opciones.tipo,'Zonal')
    Opciones.tipo='Z';
    fprintf('         > Zonal section\n');
elseif strcmp(Opciones.tipo,'zonal')
    Opciones.tipo='Z';
    fprintf('         > Zonal section\n');
elseif strcmp(Opciones.tipo,'z')
    Opciones.tipo='Z';
    fprintf('         > Zonal section\n');
elseif strcmp(Opciones.tipo,'Meridional')
    Opciones.tipo='M';
    fprintf('         > Meridional section\n');
elseif strcmp(Opciones.tipo,'meridional')
    Opciones.tipo='M';
    fprintf('         > Meridional section\n');
elseif strcmp(Opciones.tipo,'m')
    Opciones.tipo='M';
    fprintf('         > Meridional section\n');
else
    Opciones.tipo='Z';
    fprintf('         > Zonal section [Default]\n');
end
if ~isfield(Opciones,'tipo')
    Opciones.tipo='Z';
    fprintf('         > Zonal section [Default]\n');
end


if ~isfield(Opciones,'xrange')||length(Opciones.xrange)==1
    Opciones.xrange=extrem(xvar);
end
if ~isfield(Opciones,'yrange1')||length(Opciones.yrange1)==1||length(Opciones.yrange2)==1
    Opciones.yrange1=[-400 0];
end
if ~isfield(Opciones,'yrange2')||length(Opciones.yrange2)==1
    Opciones.yrange2=[nanmin(yvar(:)) -400];
end
if isfield(Opciones,'ziso')&&length(Opciones.ziso)==1
    Opciones.ziso=nice(extrem(zvar),1);
end
if ~isfield(Opciones,'zrange')||length(Opciones.zrange)==1
    Opciones.zrange=nice(extrem(zvar),12);
    if isfield(Opciones,'ziso')
        Opciones.zrange=sort([Opciones.zrange Opciones.ziso]);
    end
end
if ~isfield(Opciones,'titulo')
    Opciones.titulo='';
end
if ~isfield(Opciones,'xtitulo')
    Opciones.xtitulo='';
end
if ~isfield(Opciones,'ytitulo')
    Opciones.ytitulo='';
end
if ~isfield(Opciones,'ztitulo')
    Opciones.ztitulo='';
end
if ~isfield(Opciones,'stations')
    Opciones.stations='';
end
if ~isfield(Opciones,'colormap')
    Opciones.colormap='jet';
end

colormap(Opciones.colormap)

%% Superficie
subplot('position',[0.13 0.74 0.65 0.20]);ax1=gca;
contourf(xvar,yvar,zvar,Opciones.zrange,'LineStyle','none','clip','on');hold on;

set(ax1,'Xtick',[]);

if isfield(Opciones,'ziso')
    if ~isempty(Opciones.ziso)
        [c,h]=contour(xvar,yvar,zvar,Opciones.ziso,'k');
        clabel(c,Opciones.ziso,'fontsize',10,'color','k','rotation',0,'BackgroundColor','w','margin',1);
    end
end
if isfield(Opciones,'ziso2')
    if ~isempty(Opciones.ziso2)
        [c2,h2]=contour(xvar,yvar,zvar,Opciones.ziso2,'-','color',[0.45 0.45 0.45], ...
            'linewidth',2,'clip','on');
        clabel(c2,Opciones.ziso2,'fontsize',10,'color',[0.45 0.45 0.45], ...
            'rotation',0,'BackgroundColor','w','margin',1);
    end
end
axis([Opciones.xrange Opciones.yrange1]);
plot(xvar,ones(size(xvar))*min(Opciones.yrange1),'ko','MarkerSize',3,'MarkerEdgeColor','k', ...
    'MarkerFaceColor','k')
%Anado etiquetas a algunas estaciones
if ~isempty(Opciones.stations)
    for ilabelstations=1:length(Opciones.stations)
        text(xvar(ilabelstations),min(Opciones.yrange1),num2str(Opciones.stations(ilabelstations)), ...
            'VerticalAlignment','top','HorizontalAlignment','center','FontSize',12);
    end
end
caxis(extrem(Opciones.zrange))
%Anado batimetria
if isfield(Opciones,'batymetry')
    if exist(Opciones.batymetry,'file')>0
        B=load(Opciones.batymetry);
        if strcmp(Opciones.tipo,'Z')
            area(ax1,B.batylon-360,B.elevations,min(Opciones.yrange1),'facecolor','k','clip','on');
        else
            area(ax1,B.batylat,B.elevations,min(Opciones.yrange1),'facecolor','k','clip','on');
        end
    else
    end
end
title(sprintf('%s %s',Opciones.titulo,Opciones.ztitulo),'FontSize',12,'Fontweight','bold','interpreter','none');
ylabel(Opciones.ytitulo,'FontSize',12);

% Profundo
subplot('position',[0.13 0.1 0.65 0.57]);ax2=gca;
contourf(xvar,yvar,zvar,Opciones.zrange,'LineStyle','none');hold on;
set(ax1,'Xtick',[]);
if isfield(Opciones,'ziso')
    if ~isempty(Opciones.ziso)
        [c,h]=contour(xvar,yvar,zvar,Opciones.ziso,'k');
        clabel(c,Opciones.ziso,'fontsize',10,'color','k','rotation',0,'BackgroundColor','w','margin',1);
    end
end
if isfield(Opciones,'ziso2')
    if ~isempty(Opciones.ziso)
        [c2,h2]=contour(xvar,yvar,zvar,Opciones.ziso2,'-','color',[0.45 0.45 0.45],'linewidth',2);
        clabel(c2,Opciones.ziso2,'fontsize',10,'color',[0.45 0.45 0.45],'rotation',0,...
            'BackgroundColor','w','margin',1);
    end
end
axis([Opciones.xrange Opciones.yrange2]);

%Add estaciones
plot(xvar,ones(size(xvar))*max(Opciones.yrange2),'ko','MarkerSize',3,'MarkerEdgeColor','k',...
    'MarkerFaceColor','k')

caxis(extrem(Opciones.zrange))

%Anado batimetria
if isfield(Opciones,'batymetry')
    if exist(Opciones.batymetry,'file')>0
        B=load(Opciones.batymetry);
        if strcmp(Opciones.tipo,'Z')
            area(ax2,B.batylon-360,B.elevations,-5500,'facecolor','k');
        else
            area(ax2,B.batylat,B.elevations,min(Opciones.yrange2),'facecolor','k','clip','on');
        end
    else
    end
end

ylabel(Opciones.ytitulo,'FontSize',12)
xlabel(Opciones.xtitulo,'FontSize',12)

ax3=colorbar;
set(ax1,'position',[0.13 0.715  0.67  0.20],'FontSize',12);
set(ax2,'position',[0.13 0.100  0.67  0.58],'FontSize',12);
set(ax3,'position',[0.83 0.100  0.06  0.82],'FontSize',12);

fprintf('     >%s %s, Max:%6.3f, Min:%6.3f\n',Opciones.titulo,Opciones.ztitulo,max(zvar(:)),min(zvar(:)))

if strcmp(Opciones.tipo,'T')
    datetick('keeplimits')
end

if exist('./Secciones','dir')==0
    mkdir('.Secciones')
end

orient landscape;
CreaFigura(gcf,deblankT(fullfile('Secciones',strcat(Opciones.titulo,Opciones.ztitulo))),Opciones.figuras)
end

%%
function x = nice(v,fac)
% NICE  'nice' scaling vector.
%  X = NICE(V), where V is a 2-element vector [A,B], returns
%      an arbitrarily-but-nicely spaced vector X such that [A,B]
%      is surrounded by the limits of X.
%  X = NICE(V,F) , where F is an integer, increases the resolution
%      of X by a factor F.
%
%      V may also be a NxM matrix; in this case, only its minimum
%      and maximum value are used.
%      If minimum and maximum of V are equal, a warning message is
%      given.
%  Example:
%	           x = nice( [0.72 2.39] )
%	           y = nice( peaks )
%  returns
%	           x == [ 0.6  0.8  1.0  1.2  1.4  1.6  1.8  2.0  2.2  2.4 ] ,
%	           y == [ -8  -6  -4  -2  0  2  4  6  8  10 ] .
%  (c) 08-02-1995  Wolfgang Erasmi, Dept. Marine Physics, IfM Kiel, Germany
%      werasmi@ifm.uni-kiel.d400.de
%	Revised  31-08-1995  W. Erasmi
%           Changes to make function work even for max == min
%           and size(V) == [1 1].

% constants
mlbd1 = 4;		% boundary for 1st increase of no.s of levels
mfac1 = 2;		%     increasing factor for this case
mlbd2 = 1.8;	% boundary for 2nd increase of no.s of levels
mfac2 = 5;		%     increasing factor for this case

% error handling & defaults
if nargin < 2
    mult = 1;		% default multiplier
else
    mult = nanmax(fix(fac),1);
end

a = min(min(v));    			% minimum
b = max(max(v));  			% maximum

if ( a == b )
    %	fprintf(1,'%c',7);	% beep
    disp(' ');
    disp('nice: Warning: Upper equals lower boundary.');
    if a == 0
        a = -1;	b = 1;
    else
        a = a - abs(a);
        b = b + abs(b);
    end
end

d = b-a;% difference
od = round(log10(d)-0.5);	% order of magnitude of d; this looks somewhat
% weird but has to be used instead of fix(log10(d))
% because of numerical inaccuracy in log10.
scal = 10^od;
ds = d/scal; % relative delta (between 0.1 and 1)

% calculate stepwidth of vector; use stepwise increase of resolution,
%                                if necessary
if ds > mlbd1
    dlev = 1/mult * scal;
elseif ds > mlbd2
    dlev = 1/(mult*mfac1) * scal;
else
    dlev = 1/(mult*mfac2) * scal;
end
% create output vector
lo = floor(a/dlev)*dlev;
hi =  ceil(b/dlev)*dlev;
x = (lo:dlev:hi);
end

%%
function [X] = extrem(P)
mi = nanmin(P(:));
ma = nanmax(P(:));
X = [mi ma];
end

%%
function s=deblankT(x)
s=x(isspace(x)==0);
end
