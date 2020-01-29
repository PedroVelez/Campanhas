
%linspace(x,x0,10)

%Calculos en longitud
x= -12.4678; y=28.5227;
x0=-13.7051;y0=28.8024;
m=(y-y0)/(x-x0);
x0=-13.8933;y0=28.2256;

xf=[-13.0563 -13.1629 -13.2695 -13.3752 -13.4690 -13.5822 -13.6735 -13.7674 -13.8349 -13.8933];
yf=xf.*m-x0*m+y0;
for ii=1:length(yf)
       fprintf('%03d       ,%7.4f,   %7.4f,1\n ',24+ii,xf(ii),yf(ii))
end

%Calculos en latitud
x0=-14.0000;y0=26.5000;
x= -15.5706; y=27.7053;

m=(y-y0)/(x-x0);



yf=[ 27.7053 27.6090 27.4904 27.3393 27.2276 27.0783 26.9476 26.8724 26.8009 26.7293 26.6691 26.5861];
yf=fliplr(yf)
xf=x0+(yf-y0)/m;

for ii=1:length(yf)
       fprintf('%03d       ,%7.4f,   %7.4f,1\n ',52+ii,xf(ii),yf(ii))
end

plot(x,y,'sr');hold on
plot(x0,y0,'or');
plot(xf,yf,'ob-');

