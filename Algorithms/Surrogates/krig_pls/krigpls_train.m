function [KrigPLSmodel,Options] = krigpls_train(xtr, ytr, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coded and assembled by 
%   Ahsanul Habib
%   ahsanul.habib@student.adfa.edu.au
%   Last update August 17, 2016
%
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

%% Check for other user inputs for option
i = 1;
while i <= length(varargin)
    if strcmpi(varargin{i},'FIT_Fn') || strcmpi(varargin{i},'fitfunc') || strcmpi(varargin{i},'fit_func')
        
        FIT_Fn = varargin{i+1};
        i = i+2;
    elseif strcmpi(varargin{i},'KRG_RegressionModel') || strcmpi(varargin{i},'RegressionModel')...
            || strcmpi(varargin{i},'Regression_Model') || strcmpi(varargin{i},'Regr_Model') || strcmpi(varargin{i},'Reg_Model')...
            || strcmpi(varargin{i},'RegModel') || strcmpi(varargin{i},'RegrModel')
        
        RegressionModel = varargin{i+1};
        i = i+2;
    elseif strcmpi(varargin{i},'KRG_CorrelationModel') || strcmpi(varargin{i},'CorrelationModel')...
            || strcmpi(varargin{i},'Correlation_Model') || strcmpi(varargin{i},'Corr_Model') || strcmpi(varargin{i},'Cor_Model')...
            || strcmpi(varargin{i},'CorModel') || strcmpi(varargin{i},'CorrModel')
        
        CorrelationModel = varargin{i+1};
        i = i+2;
    elseif strcmpi(varargin{i},'KRG_Theta0') || strcmpi(varargin{i},'Theta0')...
            || strcmpi(varargin{i},'Theta_0') || strcmpi(varargin{i},'Theta_init') || strcmpi(varargin{i},'Init_theta')
        
        Theta0 = varargin{i+1};
        i = i+2;
    elseif strcmpi(varargin{i},'KRG_LowerBound') || strcmpi(varargin{i},'ThetaLB') || strcmpi(varargin{i},'Theta_LB')...
            || strcmpi(varargin{i},'Theta_LowerBound') || strcmpi(varargin{i},'Theta_Lower_Bound')
        
        Theta_LB = varargin{i+1};
        i = i+2;
    elseif strcmpi(varargin{i},'KRG_UpperBound') || strcmpi(varargin{i},'ThetaUB') || strcmpi(varargin{i},'Theta_UB')...
            || strcmpi(varargin{i},'Theta_UpperBound') || strcmpi(varargin{i},'Theta_Upper_Bound')
        
        Theta_UB = varargin{i+1};
        i = i+2;
    elseif strcmpi(varargin{i},'KRG_NumPCs') || strcmpi(varargin{i},'A') || strcmpi(varargin{i},'PCs') || strcmpi(varargin{i},'PrincipalComponents') || ...
            strcmpi(varargin{i},'Components') ||  strcmpi(varargin{i},'PrinComps') || strcmpi(varargin{i},'PComps') || ...
            strcmpi(varargin{i},'NumPCs') || strcmpi(varargin{i},'NumPrincipalComponents') || ...
            strcmpi(varargin{i},'NumComponents') || strcmpi(varargin{i},'NumPrinComps') || ...
            strcmpi(varargin{i},'NumPComps') || strcmpi(varargin{i},'Num_Components') || ...
            strcmpi(varargin{i},'Num_Princomps') || strcmpi(varargin{i},'Num_PrincipalComponents')|| ...
            strcmpi(varargin{i},'Num_Principal_Components') || strcmpi(varargin{i},'Num_PComps') || ...
            strcmpi(varargin{i},'Num_PCs') || strcmpi(varargin{i},'Max_PCs') || strcmpi(varargin{i},'MaxPCs') || ...
            strcmpi(varargin{i},'Max_PrincipalComponents') || strcmpi(varargin{i},'MaxPrincipalComponents') || ...
            strcmpi(varargin{i},'Max_Components') || strcmpi(varargin{i},'MaxComponents') || ...
            strcmpi(varargin{i},'Max_PrinComps') || strcmpi(varargin{i},'MaxPrinComps') || ...
            strcmpi(varargin{i},'Max_PComps') || strcmpi(varargin{i},'MaxPComps') || strcmpi(varargin{i},'NumComps') || strcmpi(varargin{i},'NComps')
            
        NumPCs = varargin{i+1};
        [~,~,~,~,~,~,~,Stats] = plsregress(xtr, ytr, NumPCs);
        PLS_Stats = Stats;
        i = i+2;
    else
        i = i+1;
    end
end

if ~exist('FIT_Fn','var')
    FIT_Fn = [];
end
if ~exist('RegressionModel','var')
    RegressionModel = [];   
end
if ~exist('CorrelationModel','var')
    CorrelationModel = [];
end
if ~exist('Theta0','var')
    Theta0 = [];
end
if ~exist('Theta_LB','var')
    Theta_LB = [];
end
if ~exist('Theta_UB','var')
    Theta_UB = [];
end
if ~exist('NumPCs','var')
    [~,~,~,~,~,PCTVAR,~,~] = plsregress(xtr, ytr);
    CSum = cumsum(PCTVAR(2,:)/sum(PCTVAR(2,:)));
    NumPCs = find(CSum>=0.9999 == 1, 1);
    [XL,YL,XS,YS,Beta,~,~,Stats] = plsregress(xtr, ytr, NumPCs);
    ypred_pls = [ones(size(xtr,1),1) xtr]*Beta;
    xfit = XS(:,1:NumPCs);
    yfit = Stats.Yresiduals;
    PLS_Stats = Stats;
    PLS_Stats.W = PLS_Stats.W(:,1:NumPCs);
    PLS_Stats.Xl = XL(:,1:NumPCs);
    PLS_Stats.Yl = YL(:,1:NumPCs);
    PLS_Stats.Xs = XS(:,1:NumPCs);
    PLS_Stats.Ys = YS(:,1:NumPCs);
    PLS_Stats.Beta = Beta;
    PLS_Stats.Ypls = ypred_pls;
    PLS_Stats.Xtr = xtr;
    PLS_Stats.Ytr = ytr;
    PLS_Stats.MeanXtr = mean(xtr,1);
end

Options = KrigPLSSetOptions(xfit, yfit, FIT_Fn, RegressionModel, CorrelationModel, Theta0, Theta_LB, Theta_UB, NumPCs, PLS_Stats);
[KrigPLSmodel,State] = KrigPLSFitModel(Options);
KrigPLSmodel.perf = State.Perf;

return

function Option = KrigPLSSetOptions(X, Y, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coded and assembled by
%   Ahsanul Habib
%   ahsanul.habib@student.adfa.edu.au
%   Last update August 17, 2016
%
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
% Load options
Option.SURRO_TYPE = 'KRIG+PLS';
% Initialize options
if nargin < 2 % Return empty if either X or Y missing
    Option.X = [];
    Option.Y = [];
    
    Option.FIT_Fn = [];
    
    Option.RegressionModel  = [];
    Option.CorrelationModel = [];
    Option.Theta0           = [];
    Option.Theta_LB         = [];
    Option.Theta_UB         = [];
    Option.NumPCs           = [];
    Option.PLS_Stats        = [];
    
else          % Load default values
    [npoints,nvariables] = size(X);
    
    Option.X = X;
    Option.Y = Y;
    
    Option.FIT_Fn = @krigpls_fit;
    
    Option.RegressionModel  = @krigpls_regpoly0;
    Option.CorrelationModel = @krigpls_corrgauss;
    Option.Theta0           = (npoints^(-1/nvariables))*ones(1, nvariables);
    Option.Theta_LB         = 1e-5*ones(1, nvariables);
    Option.Theta_UB         = 100*ones(1, nvariables);
    Option.NumPCs           = min(npoints-1,nvariables);
    Option.PLS_Stats        = [];
end

for i = 1:numel(varargin)
    if ~isempty(varargin{i})
        if i == 1
            Option.FIT_Fn = varargin{i};
        elseif i == 2
            Option.RegressionModel = varargin{i};
        elseif i == 3
            Option.CorrelationModel = varargin{i};
        elseif i == 4
            Option.Theta0 = varargin{i};
        elseif i == 5
            Option.Theta_LB = varargin{i};
        elseif i == 6
            Option.Theta_UB = varargin{i};
        elseif i == 7
            Option.NumPCs = varargin{i};
        elseif i == 8
            Option.PLS_Stats = varargin{i};
        end
    end
end

return

function [Surrogate, State] = KrigPLSFitModel(Option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Coded and assembled by 
%   Ahsanul Habib
%   ahsanul.habib@student.adfa.edu.au
%   Last update August 17, 2016
%
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

switch func2str(Option.FIT_Fn)
    case {'dace_fit','dacefit','dace'}
        State = srgtsFitCreateState(Option);
        if isempty(Option.Theta_LB) % no optimization for theta
            [Surrogate.model, State.Perf, State.FITFuncVal] = dace_fit(...
                Option.X, ...
                Option.Y, ...
                Option.RegressionModel, ...
                Option.CorrelationModel, ...
                Option.Theta0,[],[]);
        else
            [Surrogate.model, State.Perf, State.FITFuncVal] = dace_fit(...
                Option.X, ...
                Option.Y, ...
                Option.RegressionModel, ...
                Option.CorrelationModel, ...
                Option.Theta0, ...
                Option.Theta_LB, ...
                Option.Theta_UB);
        end
        Surrogate.model.surro_type = 'dace';
    case {'krigecr_fit','krigecrfit','krigecr','krigpls_fit','krigplsfit','krigpls'}
        State = srgtsFitCreateState(Option);
        if isempty(Option.Theta_LB) % no optimization for theta
            [Surrogate.model, State.Perf, State.FITFuncVal] = krigpls_fit(...
                Option.X, ...
                Option.Y, ...
                Option.RegressionModel, ...
                Option.CorrelationModel, ...
                Option.Theta0,[],[], ...
                Option.NumPCs, ...
                Option.PLS_Stats);
        else
            [Surrogate.model, State.Perf, State.FITFuncVal] = krigpls_fit(...
                Option.X, ...
                Option.Y, ...
                Option.RegressionModel, ...
                Option.CorrelationModel, ...
                Option.Theta0, ...
                Option.Theta_LB, ...
                Option.Theta_UB, ...
                Option.NumPCs, ...
                Option.PLS_Stats);
        end
        Surrogate.model.surro_type = 'krigpls';
end

return

function srgtSTT = srgtsFitCreateState(srgtOPT)
srgtSTT.FIT_Fn    = srgtOPT.FIT_Fn;
srgtSTT.FIT_FnVal = NaN;
return

function [Wstar,PC,alpha] = find_w_star(X,y,varargin)
%+++ Parameter settings
loop = 1;
while loop <= length(varargin)
    if strcmpi(varargin{loop},'A') || strcmpi(varargin{loop},'PCs') || strcmpi(varargin{loop},'PrincipalComponents') || ...
            strcmpi(varargin{loop},'Components') ||  strcmpi(varargin{loop},'PrinComps') || strcmpi(varargin{loop},'PComps') || ...
            strcmpi(varargin{loop},'NumPCs') || strcmpi(varargin{loop},'NumPrincipalComponents') || ...
            strcmpi(varargin{loop},'NumComponents') || strcmpi(varargin{loop},'NumPrinComps') || ...
            strcmpi(varargin{loop},'NumPComps') || strcmpi(varargin{loop},'Num_Components') || ...
            strcmpi(varargin{loop},'Num_Princomps') || strcmpi(varargin{loop},'Num_PrincipalComponents')||...
            strcmpi(varargin{loop},'num_Principal_Components') || strcmpi(varargin{loop},'Num_PComps') || ...
            strcmpi(varargin{loop},'Num_PCs')
        if ~ischar(varargin{loop+1})
            PC = varargin{loop+1};
        else
            [PC,alpha] = find_opt_pc_alpha(X,y);
        end
        loop = loop+2;
    elseif strcmp(varargin{loop},'alpha')
        if (~exist('alpha','var') || isempty(varargin{loop+1})) && ~ischar(varargin{loop+1})
            alpha = varargin{loop+1};
        end
        if exist('PC','var')
            if ~ischar(PC)
                [~,alpha] = find_opt_pc_alpha(X,y,PC);
            end
        end
        loop = loop+2;
    else
        loop = loop+1;
    end
end
if ~exist('alpha','var')
    alpha = 1;
end
[Mx,Nx] = size(X);

ssqX = sum(sum((X.^2)));  %+++ Total Variance of X
ssqY = sum(y.^2);         %+++ Total Variance of Y

T = zeros(Mx,PC);
P = zeros(Nx,PC);
W = zeros(Nx,PC);
R = zeros(1,PC);
B = zeros(Nx,PC);
R2X = zeros(1,PC);
R2Y = zeros(1,PC);

for i = 1:PC
    H = (1-alpha)*(X'*X)+alpha*(X'*y*y'*X);
    [w,eiv] = powermethod(H);
    
    t = X*w;
    p = X'*t/(t'*t);
    r = y'*t/(t'*t);
    
    W(:,i) = w;
    T(:,i) = t;
    P(:,i) = p;
    R(i) = r;
    
    B(:,i) = W(:,1:i)*(P(:,1:i)'*W(:,1:i))^(-1)*R(1:i)';
    
    R2X(i) = (T(:,i)'*T(:,i))*(P(:,i)'*P(:,i))/ssqX*100;
    R2Y(i) = (T(:,i)'*T(:,i))*(R(i)'*R(i))/ssqY*100;
    
    
    X = X-t*p';
    y = y-t*r';
end

Wstar = W*(P'*W)^(-1);  %+++ X*Wstar=T.
return

function [evec,eval] = powermethod(X)
%+++ Power method for computing inverse of a symmetric matrix.
epslon = 1e-10;
error = 1;
x0 = X(:,1);
while error > epslon
    x1 = X*x0;
    lambda = norm(x1);
    x1 = x1/lambda;
    error = norm(x1-x0)/norm(x0);
    x0 = x1;
end
evec = x1;
eval = lambda;
return

function [OptPC,OptAlpha] = find_opt_pc_alpha(X,y,PC,alpha,K,order)
%+++ Two way K-fold Cross-validation for ECR.
%+++ Input:  X: m x n  (Sample matrix)
%            y: m x 1  (measured property)
%            A: The maximal number of latent variables for cross-validation
%            K: fold. when K=m, it is leave-one-out CV
%       method: pretreatment method. Contains: center or autoscaling.
%        alpha: if = 0, ECR corresponds to PCR
%               if = 1, ECR corresponds to PLS
%          OPT: =1 : print process.
%               =0 : don't print process.
%        order: =0 : samples are ordered according to y-values.
%               =1 : samples are randomly partioned for CV.
%+++ Output: Structural data: CV

[Mx,Nx] = size(X);

if nargin < 6
    order = 0;
end
if nargin < 5
    K = 5;
end
if nargin < 4
    alpha = 1;
end
if nargin < 3
    PC = size(X,2);
end

if order == 0
    [y,index] = sort(y);
    X = X(index,:);
elseif order == 1
    X = X(randperm(size(X,1)),:);
    y = y(randperm(size(y,1)),:);
end

groups = 1+rem(0:Mx-1,K);

E = [];   %+++ store predicted values.
yytest = [];
RMSECV = zeros(PC,length(alpha));
Q2 = zeros(PC,length(alpha));
for group = 1:K
    testk = find(groups==group);
    calk = find(groups~=group);
    Xcal = X(calk,:);
    ycal = y(calk,:);
    Xtest = X(testk,:);
    ytest = y(testk,:);
    
    %+++ pretreatment
    %         [Xcal,para1,para2] = pretreat(Xcal,method);
    %         [ycal,ypara1,ypara2] = pretreat(ycal,method);
    %         Xtest = pretreat(Xtest,method,para1,para2);
    
    ypara1 = zeros(1,size(ytest,2));
    ypara2 = ones(1,size(ytest,2));
    TEMP = [];
    for a = 1:length(alpha)
        model = ECR(Xcal,ycal,PC,alpha(a));
        ypred = Xtest*model.regcoef*ypara2+ypara1;
        TEMP = [TEMP, ypred];
    end
    E = [E;TEMP];
    yytest = [yytest;ytest];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%+++ statistics for CV
for j = 1:length(alpha)
    temp = E(:,((j-1)*PC+1):(j*PC));
    error = temp - repmat(yytest,1,PC);
    PRESS = sum(error.^2);
    cvtemp = sqrt(PRESS/Mx)';
    SST = sumsqr(yytest-mean(yytest));
    
    for k = 1:PC
        SSE = sumsqr(temp(:,k)-yytest);
        Q2temp(k,1) = 1-SSE/SST;
    end
    
    RMSECV(:,j) = cvtemp;
    Q2(:,j) = Q2temp;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ecrRMSECV = min(min(RMSECV));
[ecrLV,b] = find(ecrRMSECV==RMSECV);
optAlpha = alpha(b);
OptPC = ecrLV;
OptAlpha = optAlpha;
% if ecrRMSECV <= 1e-3
%     break;
% end

return

function ecrmodel = ECR(X,y,A,alpha,method)
%+++  Elastic Component Regression, a bridge between PCR and PLS.
%+++  alpha: 0  PCR
%            1  PLS
%            (0, 1) transitional models between PCR and PLS

%+++ Parameter settings
[Mx,Nx] = size(X);
if nargin < 5
    method = 'autoscaling';
end
if nargin < 4
    alpha = 1;
end
if nargin < 3
    A = size(X,2);
end
if alpha<0 || alpha>1
    alpha = 1;
end

% [X,xpara1,xpara2] = pretreat(X,method);
% [y,ypara1,ypara2] = pretreat(y,method);
xpara1 = zeros(1,Nx);
xpara2 = ones(1,Nx);
ypara1 = zeros(1,size(y,2));
ypara2 = ones(1,size(y,2));

Xorig = X;

% A = min([Mx-1 Nx-1 A]);

ssqX = sum(sum((X.^2)));  %+++ Total Variance of X
ssqY = sum(y.^2);         %+++ Total Variance of Y

T = zeros(Mx,A);
P = zeros(Nx,A);
W = zeros(Nx,A);
R = zeros(1,A);
B = zeros(Nx,A);
R2X = zeros(1,A);
R2Y = zeros(1,A);

for i = 1:A
    H = (X'*(y*y')*X); %(1-alpha)*(X'*X)+alpha*
    [w,eiv] = powermethod(H);
    
    t = X*w;
    p = X'*t/(t'*t);
    r = y'*t/(t'*t);
    
    W(:,i) = w;
    T(:,i) = t;
    P(:,i) = p;
    R(i) = r;
    
    B(:,i) = W(:,1:i)*(P(:,1:i)'*W(:,1:i))^(-1)*R(1:i)';
    
    R2X(i) = (T(:,i)'*T(:,i))*(P(:,i)'*P(:,i))/ssqX*100;
    R2Y(i) = (T(:,i)'*T(:,i))*(R(i)'*R(i))/ssqY*100;
    
    X = X-t*p';
    y = y-t*r';
end

Wstar = W*(P'*W)^(-1);  %+++ X*Wstar=T.
% B=Wstar*R';

yfit = Xorig*B*ypara2+ypara1;

%+++ Output
ecrmodel.method = method;
ecrmodel.alpha = alpha;
ecrmodel.xpara = [xpara1;xpara2];
ecrmodel.ypara = [ypara1;ypara2];
ecrmodel.regcoef = B;
ecrmodel.yfit = yfit;
ecrmodel.wstar = Wstar;   %+++ X x wstar= T
ecrmodel.xweight = W;
ecrmodel.xscores = T;
ecrmodel.xloadings = P;
ecrmodel.yloadings = R;
ecrmodel.R2X = R2X;
ecrmodel.R2Y = R2Y;
return

% % Check input params
% if nargin == 2
%     regr = @krigpls_regpoly0;
%     corr = @krigpls_corrgauss;
%     theta0 = (m^(-1/n))*ones(1, n);
%     theta_lb = 1e-6*ones(1, n);
%     theta_ub = 200*ones(1, n);
%     num_pcs = []; ncomps = 1;
%     while true
%         [XL,YL,XS,YS,BetaVec,PCTVAR,~,STATS] = plsregress(xtr, ytr, ncomps);
%         CSum = cumsum(PCTVAR(end,:));
%         num_pcs = find(CSum>=0.95 == 1, 1);
%         if ~isempty(num_pcs)
%             break;
%         end
%         if ncomps == size(xtr, 2)
%             NumPCs = size(xtr, 2);
%             break;
%         end
%         ncomps = ncomps + 1;
%     end
%     pls_stats = STATS;
%     alphaval = 1;
% elseif nargin == 3
%     corr = @krigpls_corrgauss;
%     theta0 = (m^(-1/n))*ones(1, n);
%     theta_lb = 1e-6*ones(1, n);
%     theta_ub = 200*ones(1, n);
%     num_pcs = []; ncomps = 1;
%     while true
%         [XL,YL,XS,YS,BetaVec,PCTVAR,~,STATS] = plsregress(xtr, ytr, ncomps);
%         CSum = cumsum(PCTVAR(end,:));
%         num_pcs = find(CSum>=0.95 == 1, 1);
%         if ~isempty(num_pcs)
%             break;
%         end
%         if ncomps == size(xtr, 2)
%             NumPCs = size(xtr, 2);
%             break;
%         end
%         ncomps = ncomps + 1;
%     end
%     pls_stats = STATS;
%     alphaval = 1;
% elseif nargin == 4
%     theta0 = (m^(-1/n))*ones(1, n);
%     theta_lb = 1e-6*ones(1, n);
%     theta_ub = 200*ones(1, n);
%     num_pcs = []; ncomps = 1;
%     while true
%         [XL,YL,XS,YS,BetaVec,PCTVAR,~,STATS] = plsregress(xtr, ytr, ncomps);
%         CSum = cumsum(PCTVAR(end,:));
%         num_pcs = find(CSum>=0.95 == 1, 1);
%         if ~isempty(num_pcs)
%             break;
%         end
%         if ncomps == size(xtr, 2)
%             NumPCs = size(xtr, 2);
%             break;
%         end
%         ncomps = ncomps + 1;
%     end
%     pls_stats = STATS;
%     alphaval = 1;
% elseif nargin == 5
%     theta_lb = 1e-6*ones(1, n);
%     theta_ub = 200*ones(1, n);
%     num_pcs = []; ncomps = 1;
%     while true
%         [XL,YL,XS,YS,BetaVec,PCTVAR,~,STATS] = plsregress(xtr, ytr, ncomps);
%         CSum = cumsum(PCTVAR(end,:));
%         num_pcs = find(CSum>=0.95 == 1, 1);
%         if ~isempty(num_pcs)
%             break;
%         end
%         if ncomps == size(xtr, 2)
%             NumPCs = size(xtr, 2);
%             break;
%         end
%         ncomps = ncomps + 1;
%     end
%     pls_stats = STATS;
%     alphaval = 1;
% elseif nargin == 6
%     theta_ub = 200*ones(1, n);
%     num_pcs = []; ncomps = 1;
%     while true
%         [XL,YL,XS,YS,BetaVec,PCTVAR,~,STATS] = plsregress(xtr, ytr, ncomps);
%         CSum = cumsum(PCTVAR(end,:));
%         num_pcs = find(CSum>=0.95 == 1, 1);
%         if ~isempty(num_pcs)
%             break;
%         end
%         if ncomps == size(xtr, 2)
%             NumPCs = size(xtr, 2);
%             break;
%         end
%         ncomps = ncomps + 1;
%     end
%     pls_stats = STATS;
%     alphaval = 1;
% elseif nargin == 7
%     num_pcs = []; ncomps = 1;
%     while true
%         [XL,YL,XS,YS,BetaVec,PCTVAR,~,STATS] = plsregress(xtr, ytr, ncomps);
%         CSum = cumsum(PCTVAR(end,:));
%         num_pcs = find(CSum>=0.95 == 1, 1);
%         if ~isempty(num_pcs)
%             break;
%         end
%         if ncomps == size(xtr, 2)
%             NumPCs = size(xtr, 2);
%             break;
%         end
%         ncomps = ncomps + 1;
%     end
%     pls_stats = STATS;
%     alphaval = 1;
% elseif nargin == 8
%     [XL,YL,XS,YS,BetaVec,PCTVAR,~,STATS] = plsregress(xtr, ytr, NumComps);
%     pls_stats = STATS;
%     alphaval = 1;
% elseif nargin == 9
%     [XL,YL,XS,YS,BetaVec,PCTVAR,~,STATS] = plsregress(xtr, ytr, NumComps);
%     pls_stats = STATS;
% end