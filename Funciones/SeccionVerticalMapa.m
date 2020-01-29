function SeccionVerticalMapa(Data,Opciones)

%% Comprobaciones
fprintf('     >>>>> SeccionVertical Mapa\n');
%Data
if ~isfield(Opciones,'titulo')
	Opciones.titulo='';
end
if ~isfield(Opciones,'stations')
	Opciones.stations='';
end
if ~isfield(Opciones,'colormap')
	Opciones.colormap='jet';
end
%%
m_proj('Mercator','long',[Opciones.lon_min Opciones.lon_max],'lat',[Opciones.lat_min Opciones.lat_max]);

if ~isempty(Opciones.filecosta)
	m_usercoast(Opciones.filecosta,'patch',[.7 .6 .4,],'edgecolor',[.7 .6 .4,]);
end
for ins=1:length(Opciones.stations)
	m_line(Data.lon(ins)+360,Data.lat(ins),'marker','o','markersize',2,'color','b','MarkerFaceColor','b');
	m_text(Data.lon(ins)+360,Data.lat(ins),num2str(Opciones.stations(ins)),'Fontsize',9,'HorizontalAlignment','center','VerticalAlignment','top');
end
m_grid('linestyle','none','fontsize',08)
title(sprintf('%s',Opciones.titulo),'FontSize',12,'Fontweight','bold','interpreter','none');

if exist('./Secciones','dir')==0
	mkdir('.Secciones')
end
orient landscape;CreaFigura(gcf,deblankT(fullfile('Secciones',strcat(Opciones.titulo,Opciones.TituloSeccion,'Mapa'))),Opciones.figuras)
