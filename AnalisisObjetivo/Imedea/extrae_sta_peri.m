function extrae_sta_peri(xlon1,xlat1,NL,NC,NF,...
   xarm,yarm,zo,zint,ztot,suf,outpath);
%Esta función extrae de toda la malla de puntos, las
%estaciones de la periferia de la caja. Las guarda
%en formato *.mat por perfil, en una nueva subcarpeta
%llamada mat.
%Además, calcula todas las variables derivadas
%
%Eugenio Fraile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cargo los datos de temperatura y salinidad calculados
%por el programa de análisis objetivo
eval(['cd ',outpath,'temp\']);
load dat_temp;
eval(['cd ',outpath,'sal\']);
load dat_sal;

%genero la variable pres
pres=[zo:zint:ztot]';

%calculo las lat y lon para todos los puntos de malla
lat3=[xlat1:yarm:xlat1+((NF-1)*yarm)]';
lat2=repmat(lat3,1,NC);
lat1=reshape(lat2',1,NF*NC);

lon3=[xlon1:xarm:xlon1+((NC-1)*xarm)];
lon2=repmat(lon3,NF,1);
lon1=reshape(lon2',1,NF*NC);


%calculo el indice de la periferia de la malla
unos=ones(NF,NC);
unos(1:end,1)=nan;unos(1:end,NC)=nan;
unos(1,1:end)=nan;unos(NF,1:end)=nan;
unos=reshape(unos',1,NF*NC);
unos=isnan(unos);ind=find(unos==1);

%me quedo con las estaciones de la periferia
temp0=dat_temp(:,ind);
sal0=dat_sal(:,ind);
lat0=lat1(1,ind);
lon0=lon1(1,ind);

%grafico
%plot(lon0,lat0,'.')
%for ii=1:58;
%   text(lon0(ii),lat0(ii),num2str(ii))
%end


eval(['cd ',outpath,'mat\']);

for ii=1:size(temp0,2)
   name=[suf,sprintf('%2.2i',ii),'obj'];
   TEMP=temp0(:,ii);
   SALT=sal0(:,ii);
   LATI=lat0(1,ii);
   LONG=lon0(1,ii);
   PRES = pres;
   
   %Calculo variables derivadas
   SVAN = sw_svan(SALT,TEMP,PRES);
   xx=isnan(SVAN);
   if xx(2,:)==1
      SVAN(2,:)=SVAN(3,:);
   end
   if xx(1,:)==1
      SVAN(1,:)=SVAN(2,:);
   end
   
   dp = [PRES(1); diff(PRES)]';
   
   % make a centered specific volume anomoly
   clear csva center_p
   csva = SVAN(1,:);
   l = size(SVAN,1);
   csva(2:l) = (SVAN(1:l-1) +SVAN(2:l))/2;
   DYNH = cumsum(csva.*dp)*1e3;
   
   %make potential energy
   center_p(1) = PRES(1)/2;
   center_p(2:l) =  (PRES(1:l-1) + PRES(2:l))/2;
   g = sw_g(LATI);
   POTE = cumsum(csva.*center_p.*dp)*1e8/g;
   POTE=POTE';
   
   %make potential temp and sigma the
   PTMP = sw_ptmp(SALT,TEMP,PRES,0);
   SGTH = sw_pden(SALT,TEMP,PRES,0) - 1000;
   
   %make gamma
   LONG=360+LONG
   [GAMA,GAMA_LO,GAMA_HI] = gamma_n(SALT,TEMP,PRES,LONG,LATI);
   gama0(:,ii)=GAMA;
   
   %guardo los perfiles individualmente
   eval(['save ',name,' SALT TEMP PRES DYNH POTE PTMP SGTH GAMA LATI LONG'])
end

%ordeno N-W-S-E antes de guardar el archivo conjunto
di=diff(ind);tt=find(di>1);
tn=[size(ind,2):-1:tt(end)+2];
tw=[tt(end):-2:tt(1)];
ts=[1:tt(1)-1];
te=[tt(1)+1:2:tt(end)+1];
temp0=[temp0(:,tn) temp0(:,tw) temp0(:,ts) temp0(:,te)];
gama0=[gama0(:,tn) gama0(:,tw) gama0(:,ts) gama0(:,te)];
sal0=[sal0(:,tn) sal0(:,tw) sal0(:,ts) sal0(:,te)];
lat0=[lat0(tn) lat0(tw) lat0(ts) lat0(te)];
lon0=[lon0(tn) lon0(tw) lon0(ts) lon0(te)];

%guardo tb el archivo para toda la periferia
eval(['save ','periferia_obj',' temp0 sal0 gama0 lat0 lon0 pres'])




