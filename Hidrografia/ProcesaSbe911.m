% ProcesaSBE911 - Post SEABIRD software procesing
%
%		Este programa visualiza un perfil los perfiles del CTD sbe911
%		Despues calcula variables derivadas y lo pasa a formato .mat con
%       el mismo nombre que el fichero de entrada.
%		Es el primer paso despues de haber creado los .cnv con el software de SBE
%
%		v1.0 15 Julio 2003 - Instituto Espanol de Oceanografia
Limpia

Presmin=8;
Automatico=0;

CorreccionSalt=0; %Correccion en salinidad
EstacionesCorregirSalt=1:1:31; %Estaciones a corregir salinidad

%% Begin
if exist('../DatosCampanha.mat','file') > 0
	DC=load('../DatosCampanha.mat');
end

r=1;
np=0;
Files=dir('*.cnv');
if ~isempty(Files)
	while r==1
		np=np+1;
		keep r precorte DC Presmin Files np Automatico CorreccionSalt EstacionesCorregirSalt
		close all;
		file=Files(np).name;
		fprintf('>>>>>Procesando estacion %s  (%d/%d) \n',file,np,length(Files));
		[lati,long,nstat,gtime,depth,data,names]=Sbe911_mat(file);
		%Busco las presion
		for ivar=1:size(names,1)
			fprintf('    >Variables %s \n',names(ivar,:));
			if strfind(names(ivar,:),'Pressure, Digiquartz [db]')
				pres=data(:,ivar);
			end
		end
		presraw=pres;
		%Me quedo solo con el perfil de bajada
		iPresMax=find(data(:,1) == max(data(:,1)));
		pres=pres(1:iPresMax);
		fprintf('    >Maxima Profundidad (dbar): %d \n',max(data(:,1)));
		for ivar=1:size(names,1)
			if strfind(names(ivar,:),'Salinity, Practical [PSU]')
				saltraw=data(:,ivar);
				salt=data(1:iPresMax,ivar);
				if exist('CorreccionSalt','var')
					if ~isempty(find(EstacionesCorregirSalt==nstat, 1))
						fprintf('    >Correccion en Salinidad: %6.4f \n',CorreccionSalt)
						salt=salt+CorreccionSalt;
					end
				end
			end
			if strfind(names(ivar,:),'Salinity, Practical, 2 [PSU]')
				salt2=data(1:iPresMax,ivar);
				if exist('CorreccionSalt','var')
					if ~isempty(find(EstacionesCorregirSalt==nstat, 1))
						fprintf('    >Correccion en Salinidad2: %6.4f \n',CorreccionSalt)
						salt2=salt2+CorreccionSalt;
					end
				end
			end
			
			if strfind(names(ivar,:),'Temperature [ITS-90, deg C]')
				tempraw=data(:,ivar);
				temp=data(1:iPresMax,ivar);
			end
			if strfind(names(ivar,:),'Temperature, 2 [ITS-90, deg C]')
				temp2=data(1:iPresMax,ivar);
			end
			if strfind(names(ivar,:),'Fluorescence, WET Labs ECO-AFL/FL')
				flu=data(1:iPresMax,ivar);
			end
			if strfind(names(ivar,:),'Oxygen, SBE 43 [ml/l]')
				oxyi=data(1:iPresMax,ivar);
			end
		end
		
		%Interpolo a todos los niveles.
		temp=interp1(pres,temp,(Presmin:2:max(pres))');
		salt=interp1(pres,salt,(Presmin:2:max(pres))');
		if exist('temp2','var')
			temp2=interp1(pres,temp2,(Presmin:2:max(pres))');
			salt2=interp1(pres,salt2,(Presmin:2:max(pres))');
		end
		if exist('oxyi','var')
			oxyi=interp1(pres,oxyi,(Presmin:2:max(pres))');
		end
		pres=(Presmin:2:max(pres))';
		
		%Calculo una serie de variables derivadas
		svan = sw_svan(salt,temp,pres);
		dp = [pres(1); diff(pres)]';
		% make a centered specific volume anomoly
		csva = svan(1,:);
		l = size(svan,1);
		csva(2:l) = (svan(1:l-1) +svan(2:l))/2;
		dynh = cumsum(csva.*dp)*1e3;
		dynh=dynh';
		keyboard
		%make potential energy
		center_p(1) = pres(1)/2;
		center_p(2:l) =  (pres(1:l-1) + pres(2:l))/2;
		g = sw_g(lati);
		pote = cumsum(csva.*center_p.*dp)*1e8/g;
		pote=pote';
		%make potential tempp and sigma the
		ptmp=sw_ptmp(salt,temp,pres,0);
		sgth=sw_pden(salt,temp,pres,0)-1000;
		%make gamma
		gama = gamma_n(salt,temp,pres,long,lati);
		
		%% Figuras
		if Automatico==0
			figure(1);%set(1,'position',[0.015 0.32 0.56 0.55])
			m_proj('Mercator','long',[DC.lon_min DC.lon_max],'lat',[DC.lat_min DC.lat_max]);hold on
			m_usercoast(DC.filecosta,'patch',[.5 .5 .5]);
			m_grid('box','fancy','tickdir','in');
			title(sprintf('File: %s, estacion:%3.0d, \n Lat:%5.3f, Lon:%5.3f, depth(dbar): %4.1f \n Date:%s', ...
				file,nstat,lati,long,depth,datestr(datenum(gtime(1),gtime(2),gtime(3),gtime(4),gtime(5),gtime(6)))),'interpreter','none');
			m_line(long+360,lati,'marker','o','markersize',6,'color','b','MarkerFaceColor','b');
			m_text(long-0.15+360,lati+0.1,num2str(nstat),'Fontsize',8,'HorizontalAlignment','center','VerticalAlignment','top');
			
			figure(2);%set(gcf,'position',[0.59 0.056 0.40 0.82])
			%Data plotting
			subplot(1,3,1)                            %presiones
			plot(presraw,'.-k');hold on
			plot(pres,'b.-');
			set(gca,'ydir','reverse');grid on
			title('Presion')
			subplot(1,3,2)                            %temperatura (perfil)
			plot(tempraw,presraw,'.-k');hold on;
			plot(temp,pres,'b.-');
			if exist('temp2','var')
				plot(temp2,pres,'r');
			end
			zoom on;grid on
			set(gca,'ydir','reverse')
			title('Temperatura')
			subplot(1,3,3)                            %salinidad (perfil)
			plot(saltraw,presraw,'.-k');hold on;
			plot(salt,pres,'b.-');hold on;
			if exist('salt2','var')
				plot(salt2,pres,'r');
			end
			set(gca,'ydir','reverse');zoom on;grid on
			title('Salinidad')
			
			figure(3);%set(gcf,'position',[0.58 0.056 0.40 0.82])
			if exist('oxyi','var')
				subplot(1,3,1)
				plot(oxyi,pres)
				set(gca,'ydir','reverse');zoom on;grid on
				title('Oxygen')
				axis([2 7 -inf inf])
				
			end
			if exist('temp2','var')
				subplot(1,3,2)                            %temperatura (perfil)
				plot(temp-temp2,pres)
				set(gca,'ydir','reverse');zoom on;grid on
				title('Temperatura1-Temperatura2')
				axis([-0.002 0.002 -inf inf])
			end
			if exist('salt2','var')
				subplot(1,3,3)                            %saltinidad (perfil)
				plot(salt-salt2,pres)
				set(gca,'ydir','reverse');zoom on;grid on
				title('Salinidad1-Salinidad2')
				axis([-0.005 0.005 -inf inf])
			end
		end
		%Recorto los perfiles
		Ipm=find(pres==Presmin);
		pres=pres(Ipm:end);
		temp=temp(Ipm:end);
		salt=salt(Ipm:end);
		if exist('oxyi','var')
			oxyi=oxyi(Ipm:end);
		end
		dynh=dynh(Ipm:end);
		pote=pote(Ipm:end);
		ptmp=ptmp(Ipm:end);
		sgth=sgth(Ipm:end);
		gama=gama(Ipm:end);
		if exist('temp2','var') && exist('salt2','var')
			temp2=temp2(Ipm:end);
			salt2=salt2(Ipm:end);
		end
		%Anado las variables en formato whoi
		SALT=salt;
		TEMP=temp;
		if exist('temp2','var') && exist('salt2','var')
			TEMP2=temp2;
			SALT2=salt2;
		end
		PRES=pres;
		DYNH=dynh';
		POTE=pote;
		PTMP=ptmp;
		SGTH=sgth;
		GAMA=gama;
		if exist('oxyi','var')
			OXYI=oxyi;
		end
		LATI=lati;
		LONG=long+360;
		save(file(1:length(file)-4),'long','lati','nstat','gtime','depth','names','pres','temp','salt','sgth','ptmp','pote','gama','dynh','SALT','TEMP','PRES','DYNH','POTE','PTMP','SGTH','GAMA','LATI','LONG')
		if exist('temp2','var')
			save(file(1:length(file)-4),'-append','temp2','salt2','SALT2','TEMP2')
		end
		if exist('oxyi','var')
			save(file(1:length(file)-4),'-append','oxyi','OXYI')
		end
		fprintf('    >Datos almacenados en %s\n',strcat(file(1:length(file)-4),'.mat'))
		figure(2)
		figure(1)
		if Automatico==0 && np<length(Files)
			r=input('>>>>> Seguimos procesando ([1]/0)?');
			if r == 0
				break;
			end
			r=1;
		elseif size(Files,1)==np
			break;
		end
		clear pres oxyi dp csva center_p g pote
	end
else
	fprintf('>>>>> No hay ficheros *.cnv \n')
end
