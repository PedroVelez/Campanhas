function [x,y]=m_areatrabajo(x,y);
if nargin<2
[X,Y]=ginput(2);
[x,y]=m_xy2ll(X,Y);
end
m_line([x(1) x(2) x(2) x(1) x(1)],[y(1) y(1) y(2) y(2) y(1)],'color','k','linewidth',2)
for ii=1:2
    degrees2dms(x(ii)-360)
    degrees2dms(y(ii))
end    
return