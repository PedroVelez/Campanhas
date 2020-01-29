function make_station_grid(outpath,suf,zo,zint,ztot,NL,NC,NF);

flname=[outpath,'varobs\'];

vector=[zo:zint:ztot];
x=(ceil(NC/10))*NF;
dat1=[];
k=1;

for ii=vector
   dato=[];
   name=[flname,suf,sprintf('%4.4i',ii),'.grd'];
   fid=fopen(name,'r');
   for jj=1:5
      line=fgetl(fid);
   end
   for jj=1:x
      dat=fgetl(fid);
      dat1=str2num(dat);
      dato=merge(dato,dat1);%vector con todos los datos
   end
   dato_final(k,:)=dato;
   k=k+1;
end

var=input('** VARIABLES:   1: TEMP   2:SAL :');

if var==1
   outpath1=[outpath,'temp\'];
   dat_temp=dato_final;
   eval(['cd ',outpath1]);
   keep dat_temp
   save('dat_temp')
else
   outpath1=[outpath,'sal\'];
   dat_sal=dato_final;
   eval(['cd ',outpath1]);
   keep dat_sal
   save('dat_sal')
end
