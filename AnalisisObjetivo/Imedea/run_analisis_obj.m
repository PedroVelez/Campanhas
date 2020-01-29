function run_analisis_obj(outpath);
%Esta es la parte en la que realmente se ejecutan los
%programas fortran del IMEDEA.
%
%Eugenio Fraile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ruta=[outpath,'\programas'];
eval(['cd ',ruta]);
command = ['!makesmo'];
eval(command);

command = ['!C_pres'];
eval(command);

command = ['!analisis'];
eval(command);