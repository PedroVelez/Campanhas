Limpia


fileout='Ra1804LadcpNorte';
nstats=[1:1:24];
ANGLE=nstats.*0;

% fileout='Ra1804LadcpLaGraciosa';
% nstats=[1101:1:1110];
% ANGLE=nstats.*0;

%fileout='Ra1804LadcpCaboGhir';
%nstats=[1001:1:1016];
%ANGLE=nstats.*0;


CruiseDir='/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/2018_Raprocan1804';
LADCP_file= strcat(CruiseDir,'/LADCP/Visbeck/profiles/');

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
    file=sprintf('%sRa1804_%03d.mat',LADCP_file,nstats(ins));
    D=load(file);
    fprintf('%s\n',file)
    
    datesladcp=merge(datesladcp,datenum(D.dr.date));
    longsladcp=merge(longsladcp,D.dr.lon);
    latisladcp=merge(latisladcp,D.dr.lat);
    
    uladcp=merge(uladcp,D.dr.u);
    vladcp=merge(vladcp,D.dr.v);
    
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
    vit=interp1(p(isnan(v)==0),v(isnan(v)==0),presiladcp);
    viladcp(:,ins)=vit;
end

for ins=1:length(nstats)
    up=upladcp(:,ins);
    vp=vpladcp(:,ins);
    p=pressladcp(:,ins);
    upit=interp1(p(isnan(up)==0),up(isnan(up)==0),presiladcp);
    upiladcp(:,ins)=upit;
    vpit=interp1(p(isnan(vp)==0),vp(isnan(vp)==0),presiladcp);
    vpiladcp(:,ins)=vpit;
end

save(fileout,'longsladcp','latisladcp','upiladcp','vpiladcp','uiladcp','viladcp','uladcp','vladcp','upladcp','vpladcp','depthladcp','pressladcp','presiladcp','nstats')

