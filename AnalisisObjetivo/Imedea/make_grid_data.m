function make_grid_data...
   (outpath,xlon1,xlat1,alfo,NL,NC,NF,xarm,yarm,zo,zint,suf,studa);
%Este programa genera un archivo ascii con las variables
%introducidas en el script que necesita los programas
%fortran pora su ejecución.
%
%Eugenio Fraile
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flname=[outpath,'\info\grid.dat'];

fid=fopen(flname,'w');

fprintf(fid,'Data file about the domain and the grid, programs using this file:\n');
fprintf(fid,' C_PRES ,  ANALISIS , GEOPSUM , VAR_DERIVA , VER , OMEGA_INV .\n');
fprintf(fid,'\n');
fprintf(fid,'XLON1, XLAT1   ,  coordenates of the low left square of the domain (degrees)\n');
fprintf(fid,[num2str(xlon1),', ',num2str(xlat1),' \n']);
fprintf(fid,'\n');
fprintf(fid,'ALFO           ,  orientation: the zone can be delimited by meridians and\n');
fprintf(fid,['  ',num2str(alfo),'.           ,  parallels (ALFO=0.) or can be tilted an angle ALFO (degrees)\n']);
fprintf(fid,'\n');
fprintf(fid,' NL , NF , NC  ,  number of levels, rows and columns\n');
fprintf(fid,[' ',num2str(NL),' , ',num2str(NF),' , ',num2str(NC),'\n']);
fprintf(fid,'\n');
fprintf(fid,' XARM , YARM   ,  cell size (in degrees if ALFO=0., if not in KM\n');
fprintf(fid,[num2str(xarm),', ',num2str(yarm),'\n']);
fprintf(fid,'\n');
fprintf(fid,'Z0(STO),ZINT(STINT)    ,  upper level and level interval\n');
fprintf(fid,['   ',num2str(zo),'. ,  ',num2str(zint),'.\n']);
fprintf(fid,'\n');
fprintf(fid,' NA      STUDA       ,  root name of the used data (a4) and study area (a8)\n');
fprintf(fid,[suf,',',studa,'\n']);
fprintf(fid,'\n');
fprintf(fid,'------------------------------------------------------------------------------\n');
fprintf(fid,'\n');
fprintf(fid,'ATENTION !!!: in the chosen domain, in the graphic representation\n');
fprintf(fid,'              two files and two columns by side will be lost (the\n');
fprintf(fid,'              exterior ones):   4 files and 4 columns .\n');
fprintf(fid,'\n');
fprintf(fid,'              So it has to considereted that you must choose a bigger domain\n');
fprintf(fid,'              than the one you need.\n');
fprintf(fid,'-----------------------------------------------------------------------------\n');

fclose(fid);

eval(['cd ',outpath,'info\']);
clear fid flname
save grid 
