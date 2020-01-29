clc;close all;clear all;load Globales

Cli=load(strcat(GlobalSU.DatPath,'/Climatologias/WOA05/WOA05.mat'));

Estacion='15';

DataST{    1}=load(strcat('/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/2013_Raprocan1310/CTD/Datos/Mat/ra1310_',Estacion,'.mat'));
DataST{end+1}=load(strcat('/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/2014_Raprocan1404/CTD/Datos/Mat/ra0414_',Estacion,'.mat'));
DataST{end+1}=load(strcat('/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/2014_Raprocan1410/CTD/Datos/Mat/ra1410_',Estacion,'.mat'));
DataST{end+1}=load(strcat('/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/2015_Raprocan1502/CTD/Datos/Mat/se1502_',Estacion,'.mat'));
DataST{end+1}=load(strcat('/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/2015_Raprocan1504/CTD/Datos/Mat/ra1504_',Estacion,'.mat'));
DataST{end+1}=load(strcat('/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/2015_Raprocan1507/CTD/Datos/Mat/se1507_',Estacion,'.mat'));
DataST{end+1}=load(strcat('/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/2015_Raprocan1510/CTD/Datos/Mat/ra1510_',Estacion,'.mat'));
DataST{end+1}=load(strcat('/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/2016_Raprocan1603/CTD/Datos/Mat/ra1603_',Estacion,'.mat'));
DataST{end+1}=load(strcat('/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/2016_Raprocan1610/CTD/Datos/Mat/ra1610_',Estacion,'.mat'));
DataST{end+1}=load(strcat('/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/2017_Raprocan1704/CTD/Datos/Mat/Ra1704_',Estacion,'.mat'));
DataST{end+1}=load(strcat('/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/2017_Raprocan1710/CTD/Datos/Mat/Ra1710_',Estacion,'.mat'));

GMiL=[27.162 27.38 28.072 28.0986];
prei=1:2:4400;

Layer=[2500 3500]; %Capas donde interpolar
PTiL=[2.50 12];  %Isotermas donde iterpolar
cl=parula;

iSalCorr=[  1];
SalCorr= [0.0];

DataST{end}.salt=DataST{end}.salt2;

%% Calculos

%Datos climatologios
I=180+Locate(Cli.lon(181:end),DataST{1}.long);
J=Locate(Cli.lat,DataST{1}.lati);
CLIpre=squeeze(Cli.pre);
CLIsal=squeeze(Cli.sal(I,J,:));
CLItem=squeeze(Cli.tem(I,J,:));
CLIptm=sw_ptmp(CLIsal,CLItem,CLIpre,0);

%Aplico correcccion salinidad
DataST{iSalCorr}.salt=DataST{iSalCorr}.salt+SalCorr;

%Interpolo
for ip=1:length(DataST)
    DataST{ip}.gama2=DataST{ip}.gama.*NaN;
    ppi=DataST{ip}.pres;
    pgi=DataST{ip}.gama;
    b=diff(pgi);
    d=find(b<=0);
    ppres=ppi;
    pgama=pgi;
    disp(['d s para el perfil ',num2str(ip),': ',num2str(length(d))])
    while ~isempty(d)
        pgama(d+1)=NaN;
        ppres(d+1)=NaN;
        pptmp(d+1)=NaN;
        pgama=pgama(isnan(pgama)==0);
        ppres=ppres(isnan(ppres)==0);
        b=diff(pgama);
        d=find(b<=0);
    end
    [C,IA,IB] = intersect(ppi,ppres);
    DataST{ip}.gama2(IA)=pgama;
end

%Interpolo
for ip=1:length(DataST)
    DataST{ip}.jday=datenum(DataST{ip}.gtime);
    tempi(ip,:)=interp1(DataST{ip}.pres,DataST{ip}.temp,prei);
    salti(ip,:)=interp1(DataST{ip}.pres,DataST{ip}.salt,prei);
    presi(ip,:)=prei;
    jday(ip)=DataST{ip}.jday;
    tempiP(ip)=mean(tempi(ip,prei>Layer(1) & prei<Layer(2)));
    saltiP(ip)=mean(salti(ip,prei>Layer(1) & prei<Layer(2)));
    s=DataST{ip}.salt(:);
    ti=DataST{ip}.temp(:);
    pr=DataST{ip}.pres(:);
    Ind=find(isnan(ti)==0&isnan(s)==0&isnan(pr)==0);
    s=s(Ind);
    pt=sw_ptmp(s,ti(Ind),pr(Ind),0);
    for il=1:length(PTiL)
        saltOptmp(ip,il)=interp1(pt,s,PTiL(il));
    end
    clear s pt
end
for ip=1:1:length(DataST)
    pgama=DataST{ip}.gama2;
    pptmp=DataST{ip}.ptmp;
    psalt=DataST{ip}.salt;
    ppres=DataST{ip}.pres;
    ptmpOgama(ip,:)=interp1(pgama(isnan(pgama)==0),pptmp(isnan(pgama)==0),GMiL);
    saltOgama(ip,:)=interp1(pgama(isnan(pgama)==0),psalt(isnan(pgama)==0),GMiL);
    presOgama(ip,:)=interp1(pgama(isnan(pgama)==0),ppres(isnan(pgama)==0),GMiL);
end

color=linspace(1,64,length(jday));


%% Figuras
%% Perfiles
figure
subplot(1,2,1)
%temperatura (perfil)
for ip=1:length(DataST)
    plot(DataST{ip}.temp,DataST{ip}.pres,'.-','color',Colores(ip));hold on;grid on
end
set(gca,'ydir','reverse');zoom on;grid on
title('Temperatura')
subplot(1,2,2)
%salinidad (perfil)
for ip=1:length(DataST)
    plot(DataST{ip}.salt,DataST{ip}.pres,'.-','color',Colores(ip));hold on; grid on
end
set(gca,'ydir','reverse');zoom on;grid on
title('Salinidad')

%% Diferencias
figure
subplot(1,2,1)
for ip=2:length(DataST)
    plot(tempi(1,:)-tempi(ip,:),presi(1,:),'.-','color',Colores(ip));hold on;grid on
end
set(gca,'ydir','reverse');zoom on;grid on
title('Diferencias temperatura')

subplot(1,2,2)
for ip=2:length(DataST)
    plot(salti(1,:)-salti(ip,:),presi(1,:),'.-','color',Colores(ip));hold on;grid on
end
set(gca,'ydir','reverse');zoom on;grid on
title('Diferencias Salinidad')

%% Diagrama TS
figure
for ip=1:length(DataST)
    plot(DataST{ip}.salt,DataST{ip}.ptmp,'.-','color',cl(ceil(color(ip)),:));hold on;grid on
end
CreaFigura(strcat(mfilename,'_TS'))

%% Evolucion temporal capas
figure;
subplot(2,1,1)
plot(jday,tempiP,'-o');datetick('x',12)
title(sprintf('mean ptmp %04d-%04d',Layer(1),Layer(2)));grid on
subplot(2,1,2)
plot(jday,saltiP,'-o');datetick('x',12)
title(sprintf('mean salt %04d-%04d',Layer(1),Layer(2)));grid on
CreaFigura(strcat(mfilename,'_MeanLayer'))

figure;
for il=1:length(PTiL)
    subplot(length(PTiL),1,il)
    plot(jday,saltOptmp(:,il),'-o');hold on;datetick('x',12);grid on
    title(sprintf('salt on %6.2f isoterm',PTiL(il)))
end
CreaFigura(strcat(mfilename,'_SaltOnPT'))

figure;
for il=1:length(GMiL)
    subplot(length(GMiL),1,il)
    plot(jday,ptmpOgama(:,il),'-o');hold on;datetick('x',12);grid on
    title(sprintf('ptmp on %6.2f isopycnal',GMiL(il)))
end
CreaFigura(strcat(mfilename,'_PtmpOnGama'))

figure;
for il=1:length(GMiL)
    subplot(length(GMiL),1,il)
    plot(jday,saltOgama(:,il),'-o');hold on;datetick('x',12);grid on
    title(sprintf('salt on %6.2f isopycnal',GMiL(il)))
end
CreaFigura(strcat(mfilename,'_SaltOnGama'))
figure;
for il=1:length(GMiL)
    subplot(length(GMiL),1,il)
    plot(jday,presOgama(:,il),'-o');hold on;datetick('x',12);grid on
    title(sprintf('pres on %6.2f isopycnal',GMiL(il)))
end
CreaFigura(strcat(mfilename,'_PresOnGama'))
