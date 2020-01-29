function [x,y]=EscojeLL(n);
[X,Y]=ginput(n);
[x,y]=m_xy2ll(X,Y);
m_plot(360+x,y,'or-');
for i1=1:length(n)
    fprintf('WP%1s;%7.4f;%7.4f;1\n ',num2str(i1),x(i1)-360,y(i1))
end
return