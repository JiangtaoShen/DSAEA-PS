function  [model, perf, f] = krigpls_fit(S, Y, regr, corr, theta0, theta_lb, theta_ub, num_pcs, pls_stats)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
% Input:
% -------------------------------------------------------------------------------------------------
% S, Y          : Data points (S(i,:), Y(i,:)), i = 1,...,m
% regr          : If present, Function handle to a regression model.
%                 default - @krigecr_regpoly0
% corr          : If present, Function handle to a correlation function
%                 default - @krigecr_corrgauss.
%                 other options - @krigecr_matern_3_2 and @krigecr_matern_5_2
% theta0        : If present, Initial guess on theta, the correlation function parameters
%                 default - (m^(-1/n))*ones(1, n), where m is the no. of data
%                 points, n is the no. of variables or dimensions
% theta_lb      : If present, then lower bounds on theta
%                 default - 1e-6*ones(1, n), n is the no. of variables or
%                 dimensions. Theoretically the theta_ub is 0
% theta_ub      : If present, then upper bounds on theta
%                 default - 200*ones(1, n), n is the no. of variables or
%                 dimensions. Theoretically the theta_ub is +inf
% num_pcs       : If present, then maximum number of principal components
%                 to consider during model building. 
%
% Output:
% -------------------------------------------------------------------------------------------------
% model  : KRIGING with ECR model: a struct with the elements
%    regr       : function handle to the regression model
%    corr       : function handle to the correlation function
%    theta      : correlation function parameters
%    beta       : generalized least squares estimate
%    gamma      : correlation factors
%    sigma2     : maximum likelihood estimate of the process variance
%    S          : scaled design sites
%    Ssc        : scaling factors for design arguments
%    Ysc        : scaling factors for design ordinates
%    C          : Cholesky factor of correlation matrix
%    Ft         : Decorrelated regression matrix
%    G          : From QR factorization: Ft = Q*G'
%    OptPCs     : Optimal number of principal components
% perf  : struct with performance information. Elements
%    nv         : Number of evaluations of objective function
%    perf       : (q+2)*nv array, where q is the number of elements
%                 in theta, and the columns hold current values of
%                 [theta;  psi(theta);  type]
%               |type| = 1, 2 or 3, indicate 'start', 'explore' or 'move'
%               A negative value for type indicates an uphill step
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coded and assembled by 
%   Ahsanul Habib
%   ahsanul.habib@student.adfa.edu.au
%   Last update August 17, 2016
%% Paper Reference
% [1] MA Bouhlel, N Bartoli, A Otsmane & J Morlier (2016), "Improving kriging surrogates
%	  of high-dimensional design models by Partial Least Squares dimension reduction,"
%     Structural and Multidisciplinary Optimization, 53(5), 935-952, Chicago	
%% Code References
% [1] SN Lophaven, HB Nielsen, and J Sondergaard, "DACE - a MATLAB kriging toolbox," 
%     Tech. Rep. IMM-TR-2002-12, Technical University of Denmark, Denmark, Aug 2002, 
%     available at http://www2.imm.dtu.dk/~hbn/dace/.
% [2] FAC Viana, "SURROGATES Toolbox User's Guide," Version 3.0, 2011, 
%     available at http://sites.google.com/site/felipeacviana/surrogatestoolbox.
% [3] H Li,Q Xu and Y Liang. "libPLS: an integrated library for partial least squares
%     regression and discriminant analysis." PeerJ PrePrints 2 (2014), e190v1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Code
% Check design points
[m, n] = size(S);  % number of design sites and their dimension

sY = size(Y);
if  min(sY) == 1,
    Y = Y(:);
    lY = max(sY);
    sY = size(Y);
else
    lY = sY(1);
end

OptPCs = num_pcs;
Stats = pls_stats;

if m ~= lY
    error('S and Y must have the same number of rows');
end

% Check correlation parameters
lth = length(theta0);
if  nargin > 5  % optimization case
    if isempty(theta_lb)
        theta_lb = 1e-5*ones(1, n);
    end
    if isempty(theta_ub)
        theta_ub = 100*ones(1, n);
    end
    if  length(theta_lb) ~= lth || length(theta_ub) ~= lth
        error('theta0, lob and upb must have the same length');
    end
    if  any(theta_lb <= 0) || any(theta_ub < theta_lb)
        error('The bounds must satisfy  0 < lob <= upb');
    end
else  % given theta
    if  any(theta0 <= 0)
        error('theta0 must be strictly positive');
    end
end

% Normalize data
mS = mean(S);   sS = std(S);
mY = mean(Y);   sY = std(Y);
% 02.08.27: Check for 'missing dimension'
j = find(sS == 0);
if  ~isempty(j),  sS(j) = 1; end
j = find(sY == 0);
if  ~isempty(j),  sY(j) = 1; end
S = (S - repmat(mS,m,1)) ./ repmat(sS,m,1);
Y = (Y - repmat(mY,m,1)) ./ repmat(sY,m,1);
% [Wstar,OptPCs,OptAlpha] = find_w_star(S,Y,'num_pcs',num_pcs,'alpha',alphaval);
theta0 = theta0(1:OptPCs);
theta_lb = theta_lb(1:OptPCs);
theta_ub = theta_ub(1:OptPCs);
% Calculate distances D between points
idvec = repmat(1:size(S,1),size(S,1),1);
idmat1 = nonzeros(tril(idvec,-1));
idmat2 = nonzeros(triu(idvec,1)');
ij = [idmat1,idmat2];
D = S(idmat1,:) - S(idmat2,:);
if  min(sum(abs(D),2) ) == 0
    SS = S; YY = Y;
    % warning('Multiple design sites are not allowed. Selecting unique design sites...'),
    [uS,id_S] = unique(SS,'rows','stable');
    uY = YY(id_S,:);
    [~,dist] = knnsearch(uS,uS,'k',2);
    % nnID = id_S(dist(:,2) >= 1e-5);
    S = uS(dist(:,2) >= 1e-5,:);
    Y = uY(dist(:,2) >= 1e-5,:);
    m = size(S,1);
    % Calculate distances D between points
    idvec = repmat(1:size(S,1),size(S,1),1);
    idmat1 = nonzeros(tril(idvec,-1));
    idmat2 = nonzeros(triu(idvec,1)');
    ij = [idmat1,idmat2];
    D = S(idmat1,:) - S(idmat2,:);
end

% Regression matrix
F = feval(regr,S);  
[mF,p] = size(F);
if  mF ~= m
    error('number of rows in  F  and  S  do not match');
end
if  p > mF
    error('least squares problem is underdetermined');
end

% parameters for objective function
par = struct('corr',corr, 'regr',regr, 'y',Y, 'F',F, ...
    'D', D, 'ij',ij, 'scS',sS);

% Determine theta
if exist('theta_lb','var') || exist('theta_ub','var')
    % Bound constrained non-linear optimization
    [theta, f, fit, perf] = boxmin(theta0, theta_lb, theta_ub, par);
    if  isinf(f)
        error('Bad parameter region.  Try increasing  upb');
    end
else
    % Given theta
    theta = theta0(:);
    [f, fit] = objfunc(theta, par);
    perf = struct('perf',[theta; f; 1], 'nv',1);
    if  isinf(f)
        error('Bad point.  Try increasing theta0');
    end
end

% Return values
model = struct('Xtr',Stats.Xtr, 'Ytr',Stats.Ytr, 'regr',regr, 'corr',corr, 'theta',theta.', ...
                'beta',fit.beta, 'gamma',fit.gamma, 'sigma2',sY.^2.*fit.sigma2, ...
                'S',S, 'Ssc',[mS; sS], 'Y',Y, 'Ysc',[mY; sY], ...
                'C',fit.C, 'Ft',fit.Ft, 'G',fit.G, ...
                'OptPCs',OptPCs, 'PLSstats',Stats);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% friend functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [obj, fit] = objfunc(theta, par)
% Initialize
obj = inf;
fit = struct('sigma2',NaN, 'beta',NaN, 'gamma',NaN, ...
    'C',NaN, 'Ft',NaN, 'G',NaN);
F = par.F;
y = par.y;
D = par.D;
m = size(F,1);
% mzmx = (m^2 - m) / 2;
% Set up  R
r = feval(par.corr, theta, D);
idx = find(r > 0);   o = (1 : m)';
mu = (10+m)*eps;
R = sparse([par.ij(idx,1); o], [par.ij(idx,2); o], [r(idx); ones(m,1)+mu]);
% Cholesky factorization with check for pos. def.
[C, rd] = chol(R);
if  rd,  return, end % not positive definite

% Get least squares solution
C = C';   Ft = C \ F;
[Q, G] = qr(Ft,0);
if  rcond(G) < 1e-10
    % Check   F
    if  cond(F) > 1e15
        T = 'F is too ill conditioned\nPoor combination of regression model and design sites';
        error(T);
    else  % Matrix  Ft  is too ill conditioned
        return
    end
end
Yt = C \ y;   beta = G \ (Q'*Yt);
rho = Yt - Ft*beta;  sigma2 = sum(rho.^2)/m;
detR = prod( full(diag(C)) .^ (2/m) );
obj = sum(sigma2) * detR;
if  nargout > 1
    fit = struct('sigma2',sigma2, 'beta',beta, 'gamma',rho' / C, ...
        'C',C, 'Ft',Ft, 'G',G');
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [t, f, fit, perf] = boxmin(t0, lo, up, par)
%BOXMIN  Minimize with positive box constraints

% Initialize
[t, f, fit, itpar] = start(t0, lo, up, par);
if  ~isinf(f)
    % Iterate
    p = length(t);
    if  p <= 2
        kmax = 2;
    else
        kmax = min(p,4);
    end
    for  k = 1 : kmax
        th = t;
        [t, f, fit, itpar] = explore(t, f, fit, itpar, par);
        [t, f, fit, itpar] = move(th, t, f, fit, itpar, par);
    end
end
perf = struct('nv',itpar.nv, 'perf',itpar.perf(:,1:itpar.nv));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [t, f, fit, itpar] = start(t0, lo, up, par)
% Get starting point and iteration parameters

% Initialize
t = t0(:);  lo = lo(:);   up = up(:);   p = length(t);
D = 2 .^ ((1:p)'/(p+2));
ee = find(up == lo);  % Equality constraints
if  ~isempty(ee)
    D(ee) = ones(length(ee),1);   t(ee) = up(ee);
end
ng = find(t < lo | up < t);  % Free starting values
if  ~isempty(ng)
    t(ng) = (lo(ng) .* up(ng).^7).^(1/8);  % Starting point
end
ne = find(D ~= 1);

% Check starting point and initialize performance info
[f, fit] = objfunc(t,par);   nv = 1;
itpar = struct('D',D, 'ne',ne, 'lo',lo, 'up',up, ...
    'perf',zeros(p+2,200*p), 'nv',1);
itpar.perf(:,1) = [t; f; 1];
if  isinf(f)    % Bad parameter region
    return
end

if  length(ng) > 1  % Try to improve starting guess
    d0 = 16;  d1 = 2;   q = length(ng);
    th = t;   fh = f;   jdom = ng(1);
    for  k = 1 : q
        j = ng(k);    fk = fh;  tk = th;
        DD = ones(p,1);  DD(ng) = repmat(1/d1,q,1);  DD(j) = 1/d0;
        alpha = min(log(lo(ng) ./ th(ng)) ./ log(DD(ng))) / 5;
        v = DD .^ alpha;   tk = th;
        for  rept = 1 : 4
            tt = tk .* v;
            [ff, fitt] = objfunc(tt,par);  nv = nv+1;
            itpar.perf(:,nv) = [tt; ff; 1];
            if  ff <= fk
                tk = tt;  fk = ff;
                if  ff <= f
                    t = tt;  f = ff;  fit = fitt; jdom = j;
                end
            else
                itpar.perf(end,nv) = -1;   break
            end
        end
    end % improve
    
    % Update Delta
    if  jdom > 1
        D([1 jdom]) = D([jdom 1]);
        itpar.D = D;
    end
end % free variables

itpar.nv = nv;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [t, f, fit, itpar] = explore(t, f, fit, itpar, par)
% Explore step

nv = itpar.nv;   ne = itpar.ne;
for  k = 1 : length(ne)
    j = ne(k);   tt = t;   DD = itpar.D(j);
    if  t(j) == itpar.up(j)
        atbd = 1;   tt(j) = t(j) / sqrt(DD);
    elseif  t(j) == itpar.lo(j)
        atbd = 1;  tt(j) = t(j) * sqrt(DD);
    else
        atbd = 0;  tt(j) = min(itpar.up(j), t(j)*DD);
    end
    [ff, fitt] = objfunc(tt,par);  nv = nv+1;
    itpar.perf(:,nv) = [tt; ff; 2];
    if  ff < f
        t = tt;  f = ff;  fit = fitt;
    else
        itpar.perf(end,nv) = -2;
        if  ~atbd  % try decrease
            tt(j) = max(itpar.lo(j), t(j)/DD);
            [ff, fitt] = objfunc(tt,par);  nv = nv+1;
            itpar.perf(:,nv) = [tt; ff; 2];
            if  ff < f
                t = tt;  f = ff;  fit = fitt;
            else
                itpar.perf(end,nv) = -2;
            end
        end
    end
end % k

itpar.nv = nv;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [t, f, fit, itpar] = move(th, t, f, fit, itpar, par)
% Pattern move

nv = itpar.nv;   ne = itpar.ne;   p = length(t);
v = t ./ th;
if  all(v == 1)
    itpar.D = itpar.D([2:p 1]).^.2;
    return
end

% Proper move
rept = 1;
while  rept
    tt = min(itpar.up, max(itpar.lo, t .* v));
    [ff, fitt] = objfunc(tt,par);  nv = nv+1;
    itpar.perf(:,nv) = [tt; ff; 3];
    if  ff < f
        t = tt;  f = ff;  fit = fitt;
        v = v .^ 2;
    else
        itpar.perf(end,nv) = -3;
        rept = 0;
    end
    if  any(tt == itpar.lo | tt == itpar.up), rept = 0; end
end

itpar.nv = nv;
itpar.D = itpar.D([2:p 1]).^.25;

return