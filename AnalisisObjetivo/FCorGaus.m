%Esta funcion es el estimador de la funcion de correlacion verdadera a 
%partir de los datos observados.
%
%	k(1)=gamma es el coeficiente ruido señal 
%	k(2)=scl es la escala de correlacion.
%	d es la distancia


function y=FCorGasus(k,d)
y=(1/(1+k(1)))*exp((-d.*d)/(2*k(2)*k(2)));
return
