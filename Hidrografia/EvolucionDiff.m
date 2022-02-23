%Este programa hace diferentes graficas para ver la evolucion en el tiempo
%del drift en conductividad entre los sensores redundantes
clearvars ;close all;load Globales

datafile='Ra2012';

prescomp=1500; %Presion a partir de la cual se hacen las comparaciones
presavg=1500; %Presion a partir de la cual se hacen los promedios

%EstacionCambioCTD=32;
figura=[4];

%% Begin
DC=load('../DatosCampanha');

if ~exist('datafile','var')
    datafile=DC.campanhacode;
elseif isempty(datafile)
    datafile=DC.campanhacode;
end

load(datafile);

if exist('campanha','var')==0;campanha=DC.campanha;end

%% Figures
figure(1)
difsalte=salts-salt2s;
plot(difsalte,press);set(gca,'ydir','reverse');zoom on;grid on
axis([-0.005 0.005 -inf inf])
title(sprintf('Salinity differences [%s]',campanha))

figure(2)
diftempe=temps-temp2s;
plot(diftempe,press);set(gca,'ydir','reverse');zoom on;grid on
axis([-0.005 0.005 -inf inf])
title(sprintf('Temperature differences [%s]',campanha))

ind=find(press>prescomp);
ind2=find(press>presavg);
for i=1:size(temps,2)
    difsalt(i)=nanmean(difsalte(ind,i));
    stdifsalt(i)=nanstd(difsalte(ind,i));
    diftemp(i)=nanmean(diftempe(ind,i));
    stdiftemp(i)=nanstd(diftempe(ind,i));
    
    mtemps(i)=nanmean(temps(ind2,i));
    stdtemps(i)=nanstd(temps(ind2,i));
    mtemp2s(i)=nanmean(temp2s(ind2,i));
    stdtemp2s(i)=nanstd(temp2s(ind2,i));
    
    msalts(i)=nanmean(salts(ind2,i));
    stdsalts(i)=nanstd(salts(ind2,i));
    msalt2s(i)=nanmean(salt2s(ind2,i));
    stdsalt2s(i)=nanstd(salt2s(ind2,i));
end


[dates,I] = sort(dates);
nstats=nstats(I);

difsalt=difsalt(I);
diftemp=diftemp(I);
stdifsalt=stdifsalt(I);
stdiftemp=stdiftemp(I);

msalts=msalts(I);
stdsalts=stdsalts(I);
msalt2s=msalt2s(I);
stdsalt2s=stdsalt2s(I);

mtemps=mtemps(I);
stdsalts=stdsalts(I);
mtemp2s=mtemp2s(I);
stdsalt2s=stdsalt2s(I);

%% temperature
figure(3)
subplot(2,1,1)
errorbar(dates,msalts,stdsalts,'k','linewidth',1);hold on;grid on
h1s=plot(dates,msalts,'bs','markersize',8);hold on
h2s=plot(dates,msalt2s,'ro','markersize',4);
title(sprintf('Mean (%4d-bottom) Sal [%s]',presavg,campanha));xlabel('Fecha')
axis([floor(min(dates(isnan(difsalt)==0)))   ceil(max(dates(isnan(difsalt)==0))) -inf inf])
datetick('x','dd/mm','keeplimits')
for i1=1:5:length(nstats)
    if dates(i1)>floor(min(dates(isnan(difsalt)==0))) & dates(i1)<ceil(max(dates(isnan(difsalt)==0)))
        text(dates(i1),nanmin(msalts),num2str(nstats(i1)), ...
            'color','k','VerticalAlignment','top','fontsize',12);hold on
    end
end
legend([h1s h2s],'Sal','Sal2','Location','NorthWest')


subplot(2,1,2)
plot([floor(min(dates(~isnan(difsalt))))  ceil(max(dates(~isnan(difsalt))))],[0 0], ...
    '--','color',[0.65 0.65 0.65],'linewidth',2);grid on;hold on
errorbar(dates,difsalt,stdifsalt,'k','linewidth',1)
plot(dates,difsalt,'k','linewidth',2.25);
title(sprintf('Promedio (%4d-bottom) Sal-Sal2 [%s]',prescomp,campanha));xlabel('Fecha')
axis([floor(min(dates(~isnan(difsalt))))   ceil(max(dates(~isnan(difsalt)))) -inf inf])
datetick('x','dd/mm','keeplimits')
for i1=1:5:length(nstats)
    if dates(i1)>floor(min(dates(~isnan(difsalt)))) && dates(i1)<ceil(max(dates(~isnan(difsalt))))
        text(dates(i1),nanmax(difsalt),num2str(nstats(i1)),'color','k','VerticalAlignment','top','fontsize',12)
        plot(dates(i1),difsalt(i1),'ok');
    end
end

%% Salinity
figure(4)
subplot(2,1,1)
errorbar(dates,mtemp2s,stdsalts,'k','linewidth',1);hold on;grid on
h1t=plot(dates,mtemps,'bs','markersize',8);hold on
h2t=plot(dates,mtemp2s,'ro','markersize',4);

title(sprintf('Promedio (%4d-bottom) Tem [%s]',presavg,campanha));xlabel('Fecha')
axis([floor(min(dates(isnan(difsalt)==0)))   ceil(max(dates(isnan(difsalt)==0))) -inf inf])
datetick('x','dd/mm','keeplimits')
for i1=1:5:length(nstats)
    if dates(i1)>floor(min(dates(isnan(difsalt)==0))) && dates(i1)<ceil(max(dates(isnan(difsalt)==0)))
        text(dates(i1),nanmin(mtemps),num2str(nstats(i1)),'color','k','VerticalAlignment','top','fontsize',12)
        plot(dates(i1),mtemps(i1),'bo');hold on
        plot(dates(i1),mtemp2s(i1),'ro');
    end
end
legend([h1t h2t],'Tem','Tem2','Location','NorthWest')


subplot(2,1,2)
plot([floor(min(dates(isnan(diftemp)==0)))   ceil(max(dates(isnan(diftemp)==0)))],[0 0], ...
    '--','color',[0.65 0.65 0.65],'linewidth',2);grid on;hold on
errorbar(dates,diftemp,stdiftemp,'k');
plot(dates,diftemp,'k','linewidth',2.25);
title(sprintf('Promedio (%4d-bottom) Tem-Tem2 [%s]',prescomp,campanha));xlabel('Fecha')
axis([floor(min(dates(isnan(difsalt)==0)))   ceil(max(dates(isnan(difsalt)==0))) -inf inf])
datetick('x','dd/mm','keeplimits')
for i1=1:5:length(nstats)
    if dates(i1)>floor(min(dates(isnan(diftemp)==0))) & dates(i1)<ceil(max(dates(isnan(diftemp)==0)))
        text(dates(i1),nanmin(diftemp),num2str(nstats(i1)),'color','k','VerticalAlignment','top','fontsize',12)
        plot(dates(i1),diftemp(i1),'ko');
    end
end

CreaFigura(1,sprintf('EvolucionDiff_Sal-Sal2_%s',campanha),figura);
CreaFigura(2,sprintf('EvolucionDiff_Tem-Tem2_%s',campanha),figura);
CreaFigura(3,sprintf('EvolucionDiff_Mean(%4d-bottom)Sal-Sal2_%s',prescomp,campanha),figura);
CreaFigura(4,sprintf('EvolucionDiff_Mean(%4d-bottom)Tem-Tem2_%s',prescomp,campanha),figura);
