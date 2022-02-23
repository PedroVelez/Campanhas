function [x,y]=m_escojell(n);
[X,Y]=ginput(n);
[x,y]=m_xy2ll(X,Y);
m_plot(x,y,'or-');
for i1=1:n
fprintf('000;%7.4f;%7.4f;2\n',x(i1),y(i1))
end
return