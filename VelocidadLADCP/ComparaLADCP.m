Limpia

DC=load('../DatosCampanha');
load(DC.filebat)

stations=[1:1:24 1101:1:1110];
lev=10:10:1500;
iest=0;


DataDir1='/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/2018_Raprocan1804/LADCP/Visbeck/profiles/';
DataDir2='/Users/pvb/Dropbox/Oceanografia/Proyectos/ProcesadoLADCPLDEO_IX/Raprocan1804/processed/';

%% Metodo 1
fprintf('Reading data \n  ')
for i2=1:length(stations)
    station=stations(i2);
    FileData=sprintf('%sra1804_%03d.mat',DataDir1,station);
    if exist(FileData,'file')> 0
        fprintf('%2.2i, ',station)
        iest=iest+1;
        D=load(FileData);
        lons(iest)=D.dr.lon;
        lats(iest)=D.dr.lat;
        Data(iest).u=D.dr.u;
        Data(iest).v=D.dr.v;
        Data(iest).z=D.dr.z;
        Fechas(iest)=datenum(D.dr.date(1),D.dr.date(2),D.dr.date(3),D.dr.date(4),D.dr.date(5),D.dr.date(6));
        nstats(iest)=station;
        ulev(1:length(lev),iest)=interp1(Data(iest).z,Data(iest).u,lev);
        vlev(1:length(lev),iest)=interp1(Data(iest).z,Data(iest).v,lev);
    else
        fprintf('>>>> Falta el %sra1710_%03d.mat\n',DataDir,station);
    end
end

%% Metodo 2
fprintf('Reading data \n  ')
iest=0;
for i2=1:length(stations)
    station=stations(i2);
    FileData=sprintf('%s%03d.mat',DataDir2,station);
    if exist(FileData,'file')> 0
        fprintf('%2.2i, ',station)
        iest=iest+1;
        D=load(FileData);
        lons2(iest)=D.dr.lon;
        lats2(iest)=D.dr.lat;
        Data2(iest).u=D.dr.u;
        Data2(iest).v=D.dr.v;
        Data2(iest).v(abs(Data2(iest).v)>100)=NaN;
        Data2(iest).z=D.dr.z;
        Fechas2(iest)=datenum(D.dr.date(1),D.dr.date(2),D.dr.date(3),D.dr.date(4),D.dr.date(5),D.dr.date(6));
        nstats2(iest)=station;
        ulev2(1:length(lev),iest)=interp1(Data2(iest).z,Data2(iest).u,lev);
        vlev2(1:length(lev),iest)=interp1(Data2(iest).z,Data2(iest).v,lev);
    else
        fprintf('>>>> Falta el %sra1710_%03d.mat\n',DataDir,station);
    end
end


%
for i1=1:34
figure(i1)
plot(ulev(:,i1),-lev,'r','linewidth',2);hold on;grid on
plot(ulev2(:,i1),-lev,'r','linewidth',1);
plot(vlev(:,i1),-lev,'b','linewidth',2);
plot(vlev2(:,i1),-lev,'b','linewidth',1);
end

figure;
plot(nanmean((ulev-ulev2).^2,1),'r');hold on;grid on
plot(nanmean((vlev-vlev2).^2,1),'b')