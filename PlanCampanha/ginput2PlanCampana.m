function  ginput2PlanCampana(n)

 [lon,lat,button,ax]=m_ginput(n);
 
for ii=1:length(lat)
       disp(sprintf('%d;%7.4f;%7.4f;1',ii,lon(ii),lat(ii)))
end
