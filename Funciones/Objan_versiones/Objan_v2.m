function [Psi,oberror] = Objan_v2(xobs,yobs,obs,eta,corlen,xnew,ynew)
% OBJAN 2D - objective analysis (Optimal Statistical Interpolation).
%
%     Performs a linear estimation of the scalar field psi at the
%     coordinates (xnew,ynew) from observation data obs (located at
%     (xobs,yobs) ) by using minimum error variance methods.
%     A Gaussian spatial correlation with length scale corlen is assumed;
%     if corlen is a 2 element vector, different scales are used for
%     x and y direction.
%
% USAGE: [psi,error] = OBJAN(xobs,yobs,obs,eta,corlen,xnew,ynew)
%
%     Input:   xobs    - x coordinate of observation                     (v/M)
%              yobs    - y     "      "      "                           (v/M)
%              obs     - observations                                     (v/M)
%              eta     - signal to noise ratio (var(noise)/var(observations))  (s)
%              corlen  - correlation length (in units as of x and y)     (s/v)
%              xnew    - new x coordinates                               (v/M)
%              ynew    - new y coordinates                               (v/M)
%     Output:  psi     - estimated field at xnew, ynew                   (M)
%              oberror   - analysis noise-to-signal ratio (var(error of the mapped field)/var(observations)                 (M)
%
% Dimensions: (v) = vector, (s) = scalar, (M) = matrix
% Note:  OBJAN assumes a climatology = 0, i.e. mean field subtracted
%        (anomalies only).
%
% July-1995  W. Erasmi, Dept. Marine Physics, IfM Kiel  Germany
% April-2005 Pedro Velez, Instituto Español de Oceanografia, Updates and comements
% Error handling

if nargin < 7,
    error('Not enough input arguments.');
elseif any(size(xobs)~=size(yobs)) | any(size(yobs)~=size(obs)),
    error('xobs, yobs and obs must have the same size.');
elseif any(size(eta)~= [1 1]),
    error('Wrong size of eta');
elseif length(corlen(:)) > 2,
    error('corlen must have 1 or 2 elements.');
elseif (min(size(xnew)) > 1) | (min(size(ynew)) > 1),
    if any(size(xnew)) ~= any(size(ynew)),
        error('If xnew and ynew are given as matrices, they have to be the same size.');
    end;
elseif any(~finite(xobs)) | any(~finite(yobs)) | any(~finite(obs)) ...
        | any(~finite(eta)) | any(~finite(corlen)) ...
        | any(~finite(xnew)) | any(~finite(ynew)),
    error('OBJAN does not work for NaN''s or Inf''s.');
end
% Setting dummy values in case obs is empty
if isempty(obs);
    disp(' Warning: empty observation vectors.');
    xobs = 10*max(max(xnew));
    yobs = 10*max(max(ynew));
    obs = 0;
end;

%--------------------------------------------------------
%  Start
%--------------------------------------------------------
disp(sprintf('>>>>> Objan: Optimal statistical interpolation'))
disp(sprintf('     >noise to signal ratio: %5.2f, corlen:%5.2f',eta,corlen))

%Old version: [psi,rms] = objan(xobs,yobs,zobs,rmsobs,rmsclim,corlen,xnew,ynew);
%               rmsobs  is the variance of the observations
%               rmsclim is the variance of the noise

rmsobs=sqrt(eta*var(obs));
rmsclim=sqrt(var(obs));
if min(size(xnew)) == 1,    % xnew,ynew are vectors
    mm = length(xnew);
    nn = length(ynew);
    [xnew ynew] = meshgrid(xnew,ynew);  % make them matrices
else
    [nn,mm] = size(xnew);
end;

yi =obs(:);
N = length(yi);     % number of observations

x = [xnew(:) ynew(:)];
lenx = size(x,1);

sigmaK = rmsclim.^2;    % -> variance of climatology  ( ^= Rpsi(x,x) )
corlen=corlen(:);

if length(corlen)==1,
    corlen=corlen*[1;1];
end;

lk2 = corlen.^2;

% generate covariance matrix
% how depends the variation of each observation on the other obs
XOBS = xobs(:)*ones(1,N);
YOBS = yobs(:)*ones(1,N);
R = sigmaK * exp( -((XOBS-XOBS').^2/lk2(1) + (YOBS-YOBS').^2/lk2(2)) );

% generate Rn = R(noise); uncorrelated erros assumed
Rn = rmsobs.^2 * eye(N);

% make M and invert
disp(sprintf('     >Inverting a %d x %d matrix, ...',N,N))
Minv = inv(R+Rn);

% inner sum: sum( Minv(ij) * yi ) must be calculated only once
Q = yi' * Minv;

%solution loop
%text loop for estimate of time coefficient;
tic;    n=0;
for k=1:10;
    c = sigmaK * exp(-( (x(1,1) - xobs(:)).^2/lk2(1) + (x(1,2) - yobs(:)).^2/lk2(2) ));
    n = Q * c;                    % guess / sigmaK
    n = c'*Minv*c;                % -error - sigmaK
end;
ctim=toc*lenx/10;
if ctim > 180,
    disp(sprintf('     >Please wait approx %3.1f min.',ctim/60.))
else
    disp(sprintf('     >Please wait approx %3.1f sec.',ctim))
end;

c   = zeros(N,1);
Psi = zeros(lenx,1);
E   = zeros(lenx,1);
for p = 1:lenx,
    c = sigmaK * exp(-( (x(p,1) - xobs(:)).^2/lk2(1) + (x(p,2) - yobs(:)).^2/lk2(2) ));
    Psi(p) = Q * c;              % guess / sigmaK
    E(p)   = c'*Minv*c;          % -erro - sigmaK
end;

E   = sigmaK - E;

Psi = reshape(Psi,nn,mm);
oberror = reshape(E,nn,mm)/(rmsclim.^2);