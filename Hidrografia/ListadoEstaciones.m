close all;clear all

load('../DatosCampanha')

datafile='Ra1804';

load(datafile)


%Creo la matriz de presiones, igual en dimensiones a la de tems y sals
press=press*ones(1,size(temps,2));
%Y rellena de NaN donde no hay medida
press(isnan(temps))=NaN;

%la ordeno por fechas
[dates,I] = sort(dates);
longs=longs(I);
latis=latis(I);
press=press(:,I);
depths=depths(I);
nstats=nstats(I);


fid2=fopen(strcat('ListadoEstaciones',campanha,'.txt'),'w');
%fprintf(fid2,'%s \n','Estacion;           Fecha        Longitud;  Latitud;  Profundidad(dbar)');
for iest=1:size(longs')
    profes(iest)=max(press(:,iest));
    proffo(iest)=depths(iest);
    fprintf(fid2,'%d;%s;%7.4f;%7.4f;%8.2f\n',nstats(iest),datestr(dates(iest)),sign(longs(iest)).*floor(abs(longs(iest))),latis(iest),max(press(:,iest)));
    fprintf('%03d;     %s;  %7.4f;  %7.4f; %8.2f\n',nstats(iest),datestr(dates(iest)),sign(longs(iest)).*floor(abs(longs(iest))),latis(iest),max(press(:,iest)));
end
fclose(fid2);