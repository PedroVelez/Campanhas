function [corisom,coriso,dist] = Croscor(xobs,yobs,zobs,dint,distmax);
% croscor - compute spatial autocorrelations using two differetn statistical aproaches.
%
% USAGE: [corisom,coriso] = croscor(xobs,yobs,zobs,dint,distmax);
%
%       Input:
%               xobs      - Longitude
%               yobs      - Latitude
%               z      - observations
%               dint   -Intervalo de promedio en distancias
%               distmax-Lag maximo sobre el que se realiza el calculo de correlaciones.
%
%       Output:  
%               corisom - Autocovarianza multivariable normalizada es 
%                       decir considero que cada intervalo de distancia es 
%                       una variable. Y da la autovarianza normalizada de 
%                       cada tramo
%               coriso  - Es el estimador 'unbiased' de la Autocovarianza 
%                       normalizada considerando que cada intervalo tiene 
%                       la misma varianza. De este modo da un estimador 
%                       de: rho(dist) 1/(1+gamma), y por tanto se puede 
%                       calcular gamma mediante un ajuste.
%
% April-2006 Pedro Velez, Instituto Español de Oceanografia

xobs=xobs(:);
yobs=yobs(:);
z=zobs(:);

rd=pi/180;xkgla=60.*1.852; %Cte para pasar de grados a kilometros.

nind=round(distmax/dint); %Numero de intervalos
dist=0:dint:(nind-1)*dint; %vector con valores  iniciales de cada intervalo

mmax=0;
niso=zeros([1 nind]);
coriso=zeros([1 nind]);
den1=zeros([1 nind]);
den2=zeros([1 nind]);

varz=var(zobs);

for k=1:length(zobs)
    for l=1:length(zobs)
        yr=(yobs(l)-yobs(k))*xkgla;
        xr=(xobs(l)-xobs(k))*xkgla*cos(0.5*rd*(yobs(k)+yobs(l)));
        d=sqrt(xr^2+yr^2);
        m=floor(d/dint)+1;
        if m<=nind
            niso(m)=niso(m)+1;
            coriso(m)=coriso(m)+zobs(k)*zobs(l);
            den1(m)=den1(m)+zobs(k)^2;
            den2(m)=den2(m)+zobs(l)^2;
        end
    end
end

% corisom - Autocovarianza multivariable normalizada es decir considero que
% cada intervalo de distancia es una variable. Y da la autovarianza
% normalizada de cada tramo
corisom=coriso./sqrt(den1.*den2);
% coriso  - Es el estimador 'unbiased' de la Autocovarianza normalizada
% considerando que cada intervalo tiene la misma varianza.
% De este modo da un estimador de: rho(dist) 1/(1+gamma), y por
% tanto se puede calcular gamma
coriso=coriso./(niso*varz);