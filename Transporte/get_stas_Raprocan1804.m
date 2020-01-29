function [lat,lon,s,t,p,pt,sgth,gm,dh,pe] = get_stas_Raprocan1804(inpath,stations,deci)
% loads in selected  stations of ctd data (raw and derived quanitities)
% assumes up to three digit station number
%
% deci - specifies a decimation rate, e.g. 20 (dbar)
s = []; t = [] ; p = []; pt = []; sgth = []; gm = [];
dh = []; pe= [];s2 = []; t2 = [] ;
for st =1:length(stations)
    if stations(st)<900
        flname = [inpath,'Ra1804_',sprintf('%03d',stations(st)),'.mat'];
    elseif stations(st)>900 stations(st)<1000;
        flname = [inpath,'Ra1804_',sprintf('%3d',stations(st)),'.mat'];
    elseif stations(st)>1000
        flname = [inpath,'Ra1804_',sprintf('%d',stations(st)),'.mat'];
    end
    if exist(flname,'file') == 2
        fprintf('    > Loading %s \n',flname);
        load(flname);
        PRES = PRES(:); DYNH = DYNH(:); POTE = POTE(:);
        if nargin == 3
            %Assume its constant pressure spacing
            dp= diff(PRES(1:2));
            %Find appropriate step size
            ij = round(deci/dp);
            %Find starting index
            istart = find((PRES-ij) == min(abs(PRES-ij)));
            ii = istart:ij:length(SALT);
        else
            ii = ones(length(SALT),1);
        end
        s = merge(s,SALT(ii));
        t = merge(t,TEMP(ii));
        %s2 = merge(s2,SALT2(ii));
        %t2 = merge(t2,TEMP2(ii));
        p = merge(p,PRES(ii));
        %ox = merge(ox,OXYG(ii));
        if nargout > 6
            pt = merge(pt,PTMP(ii));
            sgth = merge(sgth,SGTH(ii));
            gm = merge(gm,GAMA(ii));
            dh = merge(dh,DYNH(ii));
            pe = merge(pe,POTE(ii));
            %fl = merge(fl,FLUOR(ii));
        end
        lat(st) = LATI;
        lon(st) = LONG;
    else
        fprintf('%s not found! \n',flname)
    end
end

% Check to see if all stations have same pressure vector
% Find pressure difference at top of each station
dsp = diff(p([1 2],:));

if any(dsp ~= dsp(1))
    fprintf('Warning!  Columns do not have equal pressure spacing \n')
else
    minp = min(p(~isnan(p)));
    maxp = max(p(~isnan(p)));
    pall = minp:dsp(1): maxp;
end
p=pall(1:length(t))';
%pp=ones(3000,1);
%ppp=merge(pp,p);p=ppp(:,2);

return
%resize arrays to full pressure size if needed
if length(pall) > size(t,1)
    n = length(pall)- size(t,1);
    t = [t; nan*ones(n,size(t,2))];
    s = [s; nan*ones(n,size(s,2))];
    %     ox = [ox; nan*ones(n,size(ox,2))];
    p  = [p; nan*ones(n,size(p,2))];
    if nargout > 6
        pt = [pt; nan*ones(n,size(pt,2))];
        sgth = [sgth; nan*ones(n,size(sgth,2))];
        gm = [gm; nan*ones(n,size(gm,2))];
    end
    if nargout >= 9
        dh = [dh; nan*ones(n,size(dh,2))];
        pe = [pe; nan*ones(n,size(pe,2))];
        %         fl = [fl; nan*ones(n,size(pe,2))];
    end
end

for i = 1:size(p,2);
    ok = ~isnan(p(:,i));
    if min(p(ok,i)) > minp
        %use  constant extrapolation to fill missing p's
        n = sum(pall < min(p(ok,i)));
        nt = [t(1,i)*ones(n,1); t(ok,i)];
        ns = [t(1,i)*ones(n,1); s(ok,i)];
        %nox = [t(1,i)*ones(n,1); ox(ok,i)];
        %%      nfl = [t(1,i)*ones(n,1); fl(ok,i)];
        nok = length(nt); 		%number of points in corrected column
        %t(1:nok,i) = nt(:);
        % s(1:nok,i) = ns(:);
        %ox(1:nok,i) = nox(:);
        if nargout > 6
            %npt = [t(1,i)*ones(n,1); pt(ok,i)];
            nsgth = [t(1,i)*ones(n,1); sgth(ok,i)];
            ngm = [t(1,i)*ones(n,1); gm(ok,i)];
            %pt(1:nok,i) = npt(:);
            sgth(1:nok,i) = nsgth(:);
            %gm(1:nok,i) = ngm(:);
        end
        if nargout >= 9
            % need to interpolate between 0 and minimum values for dh
            pp = [0 min(p(ok,i))];
            ndh = [tab1(pp,[0 dh(1,i)],pall(1:n))'; dh(ok,i)];
            % this is not quite correct for pot energy but should be ok to first order
            npe = [tab1(pp,[0 pe(1,i)],pall(1:n))'; pe(ok,i)];
            dh(1:nok,i) = ndh(:);
            pe(1:nok,i) = npe(:);
            %   fl(1:nok,i) = nfl(:);
        end
        fprintf('Column %d adjusted to fill shallow pressure\n',i)
    end
    
    
    if any(diff(p(ok,i)) ~= dsp(1))
        %use linear interpolation to fill missing p's
        fprintf('Column %d adjusted to standard pressure spacing\n',i)
        
        nt = tab1(p(ok,i),t(ok,i),pall,1);
        ns = tab1(p(ok,i),s(ok,i),pall,1);
        nox = tab1(p(ok,i),ox(ok,i),pall,1);
        %nfl = tab1(p(ok,i),fl(ok,i),pall,1);
        nok = length(nt); 		%number of points in corrected column
        t(1:nok,i) = nt(:);
        s(1:nok,i) = ns(:);
        %ox(1:nok,i) = nox(:);
        %fl(1:nok,i) = nfl(:);
        
        if nargout > 6
            npt = tab1(p(ok,i),pt(ok,i),pall,1);
            nsgth = tab1(p(ok,i),sgth(ok,i),pall,1);
            ngm = tab1(p(ok,i),gm(ok,i),pall,1);
            ndh = tab1(p(ok,i),dh(ok,i),pall,1);
            npe = tab1(p(ok,i),pe(ok,i),pall,1);
            
            pt(1:nok,i) = npt(:);
            sgth(1:nok,i) = nsgth(:);
            gm(1:nok,i) = ngm(:);
            dh(1:nok,i) = ndh(:);
            pe(1:nok,i) = npe(:);
        end
    end
end
p = pall;