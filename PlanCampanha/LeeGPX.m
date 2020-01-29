Limpia
P=gpxread('/Users/pvb/Dropbox/Oceanografia/Proyectos/Raprocan/Administracion/2016_Raprocan1710/PlanCampanha/12.gpx');

for ii=1:length(P.Latitude)
       fprintf('%0s;%7.4f;%7.4f;1\n ',P.Name{ii},P.Longitude(ii),P.Latitude(ii))
end