% CatWithNaNs - Script for objective data analysis and representation
%
% Pequeno script para redimensionar matrices utilizando nans
% 		Use:  B=catwithnans(A,nst,newcol)
% Input:Donde:	A	matriz original
%		newcol	nueva columna de datos
%		nst	n de orden (usualmente estaci?n) correspondiente a newcol

function B=CatWithNans(A,nst,newcol)
[mant,nant]=size(A);
m=length(newcol);
if m>mant
    addprof=NaN.*ones(m-mant,nant);
    B=cat(1,A,addprof);
else
    B=A;
end %if
flag=mant>m;
newcol=cat(1,newcol,(flag).*NaN.*ones(mant-m,1));
if nst>nant+1
    addstn=NaN.*ones((max(m,mant)),nst-nant-1);
    B=cat(2,B,addstn,newcol);
else
    B(:,nst)=newcol;
end %if
return
