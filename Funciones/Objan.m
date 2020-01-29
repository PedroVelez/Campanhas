function [Psi,oberror,obsi] = Objan(xobs,yobs,obs,gamma,scl,cutscl,xnew,ynew)

% Objan 2D - objective analysis (Optimal Statistical Interpolation).
%
%  Performs a linear estimation of the scalar field psi at the coordinates
%  (xnew,ynew) from observations obs (located at(xobs,yobs)) by using
%  Optimal Interpolation as described in Bretherton(1976) and Gomis and
%  Pedder (2005)
%
%  The optimal solution computed here is defined as that linear
%  combination of observed anomalies which minimises
%  (in a statistical sense) the deviations between the prediction and
%  the true field. Assuming that the anomaly field variance is constant
%  over the domain.
%
%  Given that the OI formulation is based on the idea that the observed
%  field can be represented as the sum of a large-scale mean field and a
%  smaller-scale anomaly field. This residual anomaly field, zobs, should
%  behave like a stationary, zero-mean random function of location
%
%  A Gaussian spatial correlation model with length scale scl is assumed;
%  The analyzed field is filtered (convoluted with the correlation
%  model), to dropout non-resolved scales included by the modelled
%  correlation, with characteristic scale cutscl.
%
% USAGE: [Psi,oberror,obsi] = Objan(xobs,yobs,obs,gamma,scl,cutscl,xnew,ynew)
%
%     Input:   xobs     - x coordinate of observations		(v)
%              yobs     - y     "      "      "				(v)
%              obs      - observed anomalies				(v)
%              gamma    - noise-to-signal ratio (var(noise)/var(Zobs))  (s)
%              scl      - correlation length (in units as of x and y)           (s)
%              cutscl   - characteristic cutoff length of the spatial filtering (s)
%              xnew     - new x coordinates                 (v/M)
%              ynew     - new y coordinates					(v/M)
%     Output:  psi      - estimated field at xnew, ynew     (M)
%              oberror  - analysis error
%                        (var(error of the mapped field)/var(Zobs)	(M)
%              obsi     - estimated field at xobs, yobs     (v)
%
% Dimensions: (v) = vector, (s) = scalar, (M) = matrix
%
% Julio-2010 Pedro Velez, Instituto Espa?ol de Oceanografia,
%                         Version 2.0,Scl, cutscl, gamma,...
% April-2005 Pedro Velez, Instituto Espa?ol de Oceanografia
%                       Updates and comements, error handling
% July-1995  W. Erasmi, Dept. Marine Physics, IfM Kiel  Germany
%
% Bretherton, F.P., Davis, R.E., Fandry, R.E., 1976. A technique for 1076
% objective analysis and design of oceanographic experiments 1077
% applied to MODE-73. Deep-Sea Res., Part 1, Oceanogr. Res. 1078
% Pap. 23, 559-582.


%preliminar checkings
if nargin < 8,
    error('Not enough input arguments.');
elseif any(size(xobs)~=size(yobs)) | any(size(yobs)~=size(obs)),
    error('xobs, yobs and obs must have the same size.');
elseif any(size(gamma)~= [1 1]),
    error('Wrong size of gamma');
elseif length(scl(:)) > 1,
    error('scl must have 1 element.');
elseif length(cutscl(:)) > 1,
    error('cutscl must have 1 element.');
elseif (min(size(xnew)) > 1) | (min(size(ynew)) > 1),
    if any(size(xnew)) ~= any(size(ynew)),
        error('If xnew and ynew are given as matrices, they have to be the same size.');
    end;
elseif any(~isfinite(xobs)) | any(~isfinite(yobs)) | any(~isfinite(obs)) ...
        | any(~isfinite(gamma)) | any(~isfinite(scl)) | any(~isfinite(cutscl)) ...
        | any(~isfinite(xnew)) | any(~isfinite(ynew)),
    error('>>>>> OBJAN does not work with NaN''s or Inf''s.');
end
if isempty(obs);
    error('>>>>> Error: empty observation vectors.');
end;
if min(size(xnew)) == 1,    % xnew,ynew are vectors
    mm = length(xnew);
    nn = length(ynew);
    [xnew ynew] = meshgrid(xnew,ynew);  % make them matrices
else
    [nn,mm] = size(xnew);
end;

%--------------------------------------------------------
%  Start
%--------------------------------------------------------
fprintf('     >>>>> Objan: Optimal statistical interpolation. Noise to signal ratio:%4.2f, Correlation scale:%4.2f, Cut-off scale:%4.2f\n',gamma,scl,cutscl);

yi =obs(:);             %Observations
x = [xnew(:) ynew(:)];  %Grid points
N = length(yi);         %number of observations
lenx = size(x,1);       %Number of grid points

fscl=cutscl/4;
scl2=scl^2;

% Generate autocorrelation matriz. Ee assume it can be obtained from
% a lag-autocorrelation homogeneous and isotropic model, approximated by
% a gauusian function of separation distances.
fprintf('         > Generting autocorrelation %d x %d matriz ...',N,N);
Xobs = xobs(:)*ones(1,N);
Yobs = yobs(:)*ones(1,N);
Dobs2=(Xobs-Xobs').^2+(Yobs-Yobs').^2;
R=exp(-0.5*Dobs2/scl2);

% Generate Rn = R(noise) with spatially uncorrelated noise assumed.
Rn = gamma * eye(N);

% Inverting the correllation matrix for the observations, including
% the noise.
fprintf(' Inverting a %d x %d matrix ...',N,N);
Minv=inv(R+Rn);

% Inner sum: sum( Minv(ij) * yi ) must be calculated only once
Q=yi'*Minv;

%In order to include a low-pass filter of the analyzed filter, the function
%used as model of the grid-stations correlation should be convoluted,
%obteing the following cts:
d1=sqrt(scl^2+fscl^2);
d2=sqrt(scl^2+2*fscl^2);
c1=(scl/d1)^2;
c2=(scl/d2)^2;

c   = zeros(N,1);
Psi = zeros(lenx,1);
E   = zeros(lenx,1);

%Solution
%First, test loop for estimate of time coefficient;
tic;n=0;
for p=1:floor(lenx/10);
    dg2=(x(p,1)-xobs(:)).^2+(x(p,2)-yobs(:)).^2;
    c =2*c1* exp(-0.5*dg2/d1^2) - c2* exp(-0.5*dg2/d2^2);
    Psi(p) = Q * c;
    E(p) = c'*Minv*c;
end;
ctim=toc*lenx/floor(lenx/10);
fprintf(' Computing %d mapped values.\n',lenx);
fprintf('         > Please wait approx %3.1f min. \n',ctim/60.);

for p = floor(lenx/10)+1:lenx,
    dg2=(x(p,1)-xobs(:)).^2+(x(p,2)-yobs(:)).^2;
    c =2*c1* exp(-0.5*dg2/d1^2) - c2* exp(-0.5*dg2/d2^2);
    Psi(p) = Q * c;              % Mapped field
    E(p)   = c'*Minv*c;          % Error
end;

E   = 1 - E;
Psi = reshape(Psi,nn,mm);
oberror = reshape(E,nn,mm);


% E=var(e)/var(observations), where e=(Ztrue-Zobs) is the error of the
% mapped field. Within the domain, var(e) may vary from one analysis point
% to another, in response to changes in the spatial distribution of
% observations in the vicinity of each new grid point, though var(e)
% can never exceed var(Zobs).

%Compute objectively anlyzed values onto the observation positions
obsi=obs.*NaN;
if nargout==3
    fprintf('         > Analysis onto the observations positions, please wait approx %3.1f min.\n',ctim/60.)
    x = [xobs(:) yobs(:)];  %Grid points
    lenx = size(x,1);       %Number of grid points
    [nn,mm] = size(xobs);
    for p = 1:lenx,
        dg2=(x(p,1)-xobs(:)).^2+(x(p,2)-yobs(:)).^2;
        c =2*c1* exp(-0.5*dg2/d1^2) - c2* exp(-0.5*dg2/d2^2);
        obsi(p) = Q * c;      % guess
    end;
    obsi = reshape(obsi,nn,mm);
end
fprintf('           Objan<<<<<\n')