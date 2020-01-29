function [x,y]=m_escojell(n);
[X,Y]=ginput(n);
[x,y]=m_xy2ll(X,Y);
m_plot(360+x,y,'or-');
%for ii=1:n
return