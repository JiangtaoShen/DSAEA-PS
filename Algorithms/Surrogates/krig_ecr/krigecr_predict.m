function [yhat, predstd] = krigecr_predict(x, Surrogate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
% Input:
% -------------------------------------------------------------------------------------------------
% x          : trial design sites with n dimensions.
%              For mx trial sites x:
%              If mx = 1, then both a row and a column vector is accepted,
%              otherwise, x must be an mx*n matrix with the sites stored
%              rowwise.
% Surrogate  : Structure with KRIGING with ECR model and model perf; see krigecr_fit
%
% Output:
% -------------------------------------------------------------------------------------------------
% yhat       : predicted response at x
% predstd    : estimated root mean squared error of the kriging predictor
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
try
    [yhat,predvar] = krigecr_predictor(x,Surrogate.model);
catch
    [yhat,predvar] = krigecr_predictor(x,Surrogate);
end
predstd = sqrt(predvar);
return

function  [yhat, predvar] = krigecr_predictor(x, model)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters
% Input:
% -------------------------------------------------------------------------------------------------
% x      : trial design sites with n dimensions.
%          For mx trial sites x:
%          If mx = 1, then both a row and a column vector is accepted,
%          otherwise, x must be an mx*n matrix with the sites stored
%          rowwise.
% model  : KRIGING with ECR model; see krigecr_fit
%
% Output:
% -------------------------------------------------------------------------------------------------
% yhat   : predicted response at x
% predvar: estimated mean squared error of the kriging predictor
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
[m, n] = size(model.S);  % number of design sites and number of dimensions
sx = size(x);            % number of trial sites and their dimension
if  min(sx) == 1 && n > 1 % Single trial point
    nx = max(sx);
    if  nx == n
        mx = 1;  
        x = x(:).';
    end
else
    mx = sx(1);  
    nx = sx(2);
end

% Normalize trial sites
x = (x - repmat(model.Ssc(1,:),mx,1)) ./ repmat(model.Ssc(2,:),mx,1);
q = size(model.Ysc,2);  % number of response functions

if  mx == 1  % one site only
    td = nan(mx*m,model.OptPCs);
    for p = 1:model.OptPCs
        dx = (repmat(x,m,1) - model.S).*repmat(model.Wstar(:,p)',m,1);  % distances to design sites
        td(:,p) = feval(model.corr, model.theta(p), dx);
    end
    f = feval(model.regr, x);
    r = prod(td, 2);
    % scaled dace_evaluate
    sy = f * model.beta + (model.gamma*r).';
    
    % prediction
    yhat = (model.Ysc(1,:) + model.Ysc(2,:) .* sy)';

    % prediction variance
    rt = model.C \ r;
    u  = model.Ft.' * rt - f.';
    v  = model.G \ u;
    
    predvar = repmat(model.sigma2,mx,1) .* repmat((1 + sum(v.^2) - sum(rt.^2))',1,q);

else  % several trial sites
    % get distances to design sites
    dx = zeros(mx*m,n);  kk = 1:m;
    for  k = 1 : mx
        dx(kk,:) = repmat(x(k,:),m,1) - model.S;
        kk = kk + m;
    end
    td = nan(mx*m,model.OptPCs);
    for p = 1:model.OptPCs
        dxx = dx.*repmat(model.Wstar(:,p)',mx*m,1);  % distances to design sites
        td(:,p) = feval(model.corr, model.theta(p), dxx);
    end
    % get regression function and correlation
    f = feval(model.regr, x);
    r = prod(td, 2);%feval(model.corr, model.theta, dx);
    r = reshape(r, m, mx);

    sy = f * model.beta + (model.gamma * r).';
    yhat = repmat(model.Ysc(1,:),mx,1) + repmat(model.Ysc(2,:),mx,1) .* sy;
    
    rt = model.C \ r;
    u = model.G \ (model.Ft.' * rt - f.');
    predvar = repmat(model.sigma2,mx,1) .* repmat((1 + colsum(u.^2) - colsum(rt.^2))',1,q);

end % of several sites

return

function  s = colsum(x)

if  size(x,1) == 1
    s = x;
else
    s = sum(x);
end

return
