%Este programa calcula transportes geostroficos a partir de los programas de WHOI
%shell to pass vars to mk_circ2.m
%
%Todos los ficheros de texto deben de ser editados acordemente para definir:
%	Directorio donde se encuentran los ficheros
%	Pares de estaciones y capas de no movimiento

output_dir='/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/2018_Raprocan1804/Analisis/Transporte/';
layer_file='sig_Raprocan1804_gamma.dat';

if 0
    %% Agadir
    %Fichero de informacion con los pares de estaciones y la capa de no-movimiento
    %Tambien se especifico el directorio donde estan los datos.
    cntrl_file='cntrl_Agadir.dat';
    %Fichero que contienen la informaci?n necesaria sobre las capas de gamma_n.
    %bot_file= 'a22.sea';
    %Fichero en el que se guardara los resultados.
    output_file=strcat(output_dir,'trans_masa_Agadir');
    %Fichero de batimetr?a
    %bath_file='bat_36N.dat';
    mk_circ2_Raprocan1804(cntrl_file,layer_file,output_file);
    
    %% LaGraciosa
    %Fichero de informaci?n con los pares de estaciones y la capa de no-movimiento
    %Tambien se especifico el directorio donde estan los datos.
    cntrl_file='cntrl_LaGraciosa.dat';
    %Fichero que contienen la informaci?n necesaria sobre las capas de gamma_n.
    %bot_file= 'a22.sea';
    %Fichero en el que se guardara los resultados.
    output_file=strcat(output_dir,'trans_masa_LaGraciosa');
    %Fichero de batimetr?a
    %bath_file='bat_36N.dat';
    mk_circ2_Raprocan1804(cntrl_file,layer_file,output_file);
    
    
    %% Norte
    %Fichero de informacion con los pares de estaciones y la capa de no-movimiento
    %Tambien se especifico el directorio donde estan los datos.
    cntrl_file='cntrl_Norte.dat';
    %Fichero que contienen la informaci?n necesaria sobre las capas de gamma_n.
    %bot_file= 'a22.sea';
    %Fichero en el que se guardara los resultados.
    output_file=strcat(output_dir,'trans_masa_Norte');
    %Fichero de batimetr?a
    %bath_file='bat_36N.dat';
    mk_circ2_Raprocan1804(cntrl_file,layer_file,output_file);
    
    %% Lanzarote
    %Fichero de informaci?n con los pares de estaciones y la capa de no-movimiento
    %Tambien se especifico el directorio donde estan los datos.
    cntrl_file='cntrl_Lanzarote.dat';
    %Fichero que contienen la informaci?n necesaria sobre las capas de gamma_n.
    %bot_file= 'a22.sea';
    %Fichero en el que se guardara los resultados.
    output_file=strcat(output_dir,'trans_masa_Lanzarote');
    %Fichero de batimetria
    %bath_file='bat_36N.dat';
    mk_circ2_Raprocan1804(cntrl_file,layer_file,output_file);
end

    %% CJuby
    %Fichero de informaci?n con los pares de estaciones y la capa de no-movimiento
    %Tambien se especifico el directorio donde estan los datos.
    cntrl_file='cntrl_CJuby.dat';
    %Fichero que contienen la informaci?n necesaria sobre las capas de gamma_n.
    %bot_file= 'a22.sea';
    %Fichero en el que se guardara los resultados.
    output_file=strcat(output_dir,'trans_masa_CJuby');
    %Fichero de batimetria
    %bath_file='bat_36N.dat';
    mk_circ2_Raprocan1804(cntrl_file,layer_file,output_file);
