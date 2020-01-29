Limpia


fileout='Ra1804Ladcp2Norte';
nstats=[1:1:4 6:1:24];
ANGLE=nstats.*0;

%fileout='Ra1710LadcpLaGraciosa';
%nstats=[1101:1:1109 1016];
%ANGLE=nstats.*0;

%fileout='Ra1710LadcpCaboGhir';
%nstats=[1001:1:1016];
%ANGLE=nstats.*0;


CruiseDir='/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/2017_Raprocan1710';
LADCP_file='/Users/pvb/Dropbox/Oceanografia/Proyectos/ProcesadoLADCPLDEO_IX/Raprocan1804/processed/';

%% Inicio
longsladcp=[];
latisladcp=[];
uladcp=[];
vladcp=[];
upladcp=[];
vpladcp=[];
depthladcp=[];
pressladcp=[];
datesladcp=[];

for ins=1:length(nstats)
    fprintf('     > Voy por la %s, ',num2str(nstats(ins)))
    file=sprintf('%s%03d.mat',LADCP_file,nstats(ins));
    D=load(file);
    fprintf('%s\n',file)
    
    datesladcp=merge(datesladcp,datenum(D.dr.date));
    longsladcp=merge(longsladcp,D.dr.lon);
    latisladcp=merge(latisladcp,D.dr.lat);
    
    up=D.dr.u;
    vp=D.dr.v;
    vp(abs(vp)>100)=NaN;
    uladcp=merge(uladcp,up);
    vladcp=merge(vladcp,vp);
    
    upladcpP=D.dr.u.*cosd(ANGLE(ins))+D.dr.v.*sind(ANGLE(ins));
    vpladcpP=D.dr.u.*sind(ANGLE(ins))+D.dr.v.*cosd(ANGLE(ins));
    
    upladcp=merge(upladcp,upladcpP);
    vpladcp=merge(vpladcp,vpladcpP);
    depthladcp=merge(depthladcp,D.dr.z);
    
    preladcp=sw_pres(D.dr.z,D.dr.lat);
    pressladcp=merge(pressladcp,preladcp);
    
end

presiladcp=8:8:ceil(max(pressladcp(:)));
for ins=1:length(nstats)
    u=uladcp(:,ins);
    v=vladcp(:,ins);
    p=pressladcp(:,ins);
    uit=interp1(p(isnan(u)==0),u(isnan(u)==0),presiladcp);
    uiladcp(:,ins)=uit;
    if length(v(isnan(v)==0))==0
        vit=presiladcp.*NaN;
    else
        vit=interp1(p(isnan(v)==0),v(isnan(v)==0),presiladcp);
    end
    viladcp(:,ins)=vit;
end

for ins=1:length(nstats)
    up=upladcp(:,ins);
    vp=vpladcp(:,ins);
    p=pressladcp(:,ins);
    
    if length(up(isnan(up)==0))==0
        upit=presiladcp.*NaN;
    else
        upit=interp1(p(isnan(up)==0),up(isnan(up)==0),presiladcp);
    end
    
    
    upiladcp(:,ins)=upit;
    
    if length(vp(isnan(vp)==0))==0
        vpit=presiladcp.*NaN;
    else
        vpit=interp1(p(isnan(vp)==0),vp(isnan(vp)==0),presiladcp);
    end
    
    
    vpiladcp(:,ins)=vpit;
end

save(fileout,'longsladcp','latisladcp','upiladcp','vpiladcp','uiladcp','viladcp','uladcp','vladcp','upladcp','vpladcp','depthladcp','pressladcp','presiladcp','nstats')

