function mk_circ2_Raprocan1804(ctrl_file,layer_file,output_file,bath_file);
%function mk_circ2(ctrl_file,layer_file,output_file,bathymetry_file);
% Computes geostrophic circulation from matlab ctd station files
% Includes transport quanitites for heat, salt and ctd oxygen
%
% LE HE QUITADO EL OXIGENO YA QUE EN LA CAMPA?A NO SE LLEVO SENSOR DE
% OXIGENO CON EL CTD
%
%  INPUTS:
%      ctrl_file: name of control files to use
%      layer_file: name of file defining layers used in analysis
%      output_file: name of file to save results to
%      bathymetry_file: name of file of bathymetric depths
%
%  FILE FORMATS
%    ctrl_file (ascii file)
%       (optional header)  all lines must begin with '%'
%       line 1: full path of location of matlab ctd station files
%       line 2: parameter used for geostrophic reference
%                 (pres, gamma, temp, ptemp, sgth)
%       line 3: number of paramters in each station entry (below)
%       line 4 and higher: one line entries for each station pair of form
%
%    station_#1 station_#2 reference_sfc triangle_method intervening_topo
%
%       triangle_method and intervening_topo are optional
%       see mk_geovel for info on triangle_method and intervening_topo
%
%
%      layer_file (ascii file)
%         (optional header)  all lines must begin with '%'%
%         line 1: parameter used to define layers
%                 (pres, gamma, temp, ptemp, sgth)
%         line 2 & higher: list of values of layer interfaces
%
%      output_file (matlab format)
%         contains variables use for analysis of circulation including the
%         total transport in each layer for each station pair and transport
%         in bottom triangles.
%
%      bathymetry_file (ascii format)
%                three column file of lat, lon and depth (m) along track
%
% copywrite Paul E. Robbins, 1995
%
% corrected 12/95 to use in situ density
% load in control file
fprintf('>>>>> %s\n',mfilename);

if nargin == 0;
    ctrl_file = input('Enter name of control file ','s');
end
fid = fopen(ctrl_file);
if fid == -1;
    disp(['Control file ',ctrl_file,' not found'])
    return
end

%load in directory where ctd stations are found
header = fgetl(fid);
while strcmp(header(1),'%')
    header = fgetl(fid);
end
indir = header; indir = deblank(indir);

%load in name of paramter to use as reference level
%should be either pres,temp,ptemp,sgth,or gamma;
ref_param  =fgetl(fid); ref_param = deblank(ref_param);
%
%read in number of controls per line
nctrl = fscanf(fid,'%i',1);
%load in block of dat defnining
%station pairs, reference levels and triangle methods;
ctrls = fscanf(fid,'%f',[nctrl,inf])';
fclose(fid);

if nargin < 2
    layer_file = input('Enter name of Layer Definition file ','s');
    layer_file = deblank(layer_file);
end
if isempty(layer_file)
    %assume just one layer
    cuts = [0 1e6];
    lay_param = 'pres';
    layer_file = 'One Layer';
else
    fid = fopen(layer_file);
    if fid == -1
        %check in the default layer_defn directory to see if its there
        fid = fopen(['/n3/layer_defn/',layer_file]);
    end
    if fid == -1;
        % if still can't find it
        fprintf('    > Layer Definition file %s not found',layer_file)
        return
    end
    header = fgetl(fid);
    while strcmp(header(1),'%')
        header = fgetl(fid);
    end
    %load in name of paramter to use define layers
    lay_param  =header; lay_param = deblank(lay_param);
    cuts = fscanf(fid,'%f',[1,inf]); cuts = sort(cuts);
    fclose(fid);
end

if nargin < 3
    output_file = input('Enter name for output file ','s');
end

nsta = size(ctrls,1);

fprintf('    > Computing transport for %d station pairs from directory %s\n',nsta,indir)
fprintf('    > Reference levels based on %s\n',ref_param)
fprintf('    > Analysis of %d layers defined by %s\n',length(cuts)-1,lay_param)
fprintf('    > Writing output to file %s \n',output_file)

if nargin == 4
    eval(['load ',bath_file])
    disp(['Using bottom bathymetry from file ',bath_file])
    %  figure out variable name of bathymetry data
    fs = find(bath_file == '/');
    if any(fs)
        bath_file = bath_file(fs(length(fs))+1:length(bath_file));
    end
    fs = find(bath_file == '.');
    if any(fs)
        bath_file = bath_file(1:fs(1)-1);
    end
    eval(['bot_topo = ',bath_file,';'])
    %convert meters to dbar
    bot_topo(:,3) = sw_pres(bot_topo(:,3),bot_topo(:,1));
end

%intialize arrays to proper size to optimize speed
area = ones(length(cuts)-1,nsta);
mass = area; tri_area = area; vol_trans = area;
tri_vol_trans = area;  mass_trans = area; tri_mass_trans = area;

% Jane's change
c1=1;c2=2;
stnvel=ones(3300,nsta)*nan;

%%
% Alonso's change. I introduce pllleea and ddeel
ddeel=10;
fprintf('    > Using a decimation rate of %d dbar',ddeel);
for s = 1:nsta
    [lat,lon,salt,temp,pres,ptemp,sgth,gamma,dh,ppllpp]=get_stas_Raprocan1804(indir,ctrls(s,[1 2]),ddeel);
    %calculate a geostrophic velocity profile
    if strcmp(ref_param,'pres')
        eval(['ref_sfc = ',ref_param,';']);
    else
        eval(['ref_sfc=transport_prop(',ref_param,');']);
    end
    %call mk_geovel with appropriate number of arguments (based on ctrl file)
    if nctrl == 3
        tri_method = 1; 		% assume constant velocity in bottom triangles
    else
        tri_method = ctrls(s,4);
    end
    if nctrl < 5
        int_topo = 9999;       % assume imposibly deep intervening topo
    else
        int_topo = ctrls(s,5);
    end
    
    if  nargin <= 3
        [stvel,starea,tri_idx,maxd(s),ref_used(s)] = ...
            mk_geovel(dh,lat,lon,ref_sfc,ctrls(s,3),pres,tri_method,int_topo);
    else
        % use bottom topography file
        sta_topo = station_topo(lat,lon,bot_topo);
        [stvel,starea,tri_idx,maxd(s),ref_used(s)] = ...
            mk_geovel(dh,lat,lon,ref_sfc,ctrls(s,3),pres,tri_method,int_topo,sta_topo);
    end
    %%%%%%%%%%%%%%%%%%%Jane's change %%%%%%%%%%%%%%%%%%
    % store station velocities
    for jj=1:length(stvel)
        nn=find(~isnan(pres));
        stnvel(nn,c1)=sw_dpth(pres(nn), mean(lat));
        stnvel(nn,c2)=stvel(nn);
    end
    c1=c1+2;
    c2=c2+2;
    %calculate transport properties
    ctemp = transport_prop(temp);
    cptemp = transport_prop(ptemp);
    crho = transport_prop(sw_dens(salt,temp,pres(:)));
    %crho = transport_prop(sgth+1000);
    %cox = transport_prop(ox).*crho.*starea/1e6;  %oxygen in moles/meter
    csalt = transport_prop(salt);
    cp = sw_cp(csalt,ctemp,pres(:));
    
    csalt = csalt.*crho.*starea/1e3; 	%salt in kg/meter
    heat = cp.*cptemp.*crho.*starea; 	%heat in J/meter
    clat(s) = mean(lat);
    % need to take into account crossing date line , added 7/17/96
    if sign(lon(1)) == sign(lon(2)) | abs(diff(lon)) < 180;
        clon(s) = mean(lon);
    else
        lon(lon<0) = lon(lon<0)+360;
        clon(s) = mean(lon);
        if clon > 180; clon = clon-360;end
    end
    eval(['lay_sfc = ',lay_param,';']);
    if strcmp(lay_param,'pres')
        ll = lay_sfc(:);
    else
        ll = transport_prop(lay_sfc);
    end
    
    for c = 1:length(cuts)-1
        fs = ll >= cuts(c) & ll < cuts(c+1);
        % now do sums of quatities for each layer
        area(c,s) = sum(starea.*fs);
        mass(c,s) = sum(starea.*fs.*crho);
        tri_area(c,s) = sum(starea.*tri_idx.*fs);
        vol_trans(c,s)=sum(starea.*stvel.*fs);
        tri_vol_trans(c,s)=sum(starea.*stvel.*tri_idx.*fs);
        mass_trans(c,s)=sum(starea.*stvel.*crho.*fs);
        tri_mass_trans(c,s)=sum(starea.*stvel.*crho.*tri_idx.*fs);
        
        %now calculate transport properties
        trans1(c,s) = sum(heat.*stvel.*fs);
        trans2(c,s) = sum(csalt.*stvel.*fs);
        %trans3(c,s) = sum(cox.*stvel.*fs);
        prop1(c,s) = sum(heat.*fs);
        prop2(c,s) = sum(csalt.*fs);
        %prop3(c,s) = sum(cox.*fs);
    end
    stlat(s) = lat(1);
    stlon(s) = lon(1);
    stnum(s) = ctrls(s,[1]);
end

stlat(s+1) = lat(2);
stlon(s+1) = lon(2);
stnum(s+1) = ctrls(s,[2]);

trans_labels = str2mat('Heat','Salinity');
trans_units = str2mat('Watts','kg/s');

eval(['save ',output_file, ' area stvel tri_area vol_trans '...
    'tri_vol_trans mass_trans tri_mass_trans ctrl_file ' ...
    'layer_file mass clat clon maxd ref_used cuts '...
    'lay_param trans_labels trans_units trans1 trans2 '...
    'prop1 prop2 stlat stlon stnum'])

velfile=[strrep(output_file,'trans_masa','vel')];
save(velfile,'stnvel')