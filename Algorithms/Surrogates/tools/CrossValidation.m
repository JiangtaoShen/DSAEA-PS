function [PRESSRMS, Rvec, model_predvar, option_predvar, eXV, predvarXVc, yhatXVc , srgtSRGTXV, srgtOPTXV] = CrossValidation(PARAM, SrgtOPT, SrgtSRGT, kfolds, idxFolds)
%Function srgtsCrossValidation computes cross-validation measures for the
%given surrogate models. Thus, for example:
%
%     PRESSRMS = srgtsCrossValidation(OPTIONS): returns a 1 x 1 vector of
%     PRESSRMS value for each surrogate, which is given by:
%
%     PRESSRMS = sqrt ( PRESS / NBPOINTS );
%
%     where:
%         * PRESS (Prediction Sum of Squares): PRESS = (eXV.')*eXV
%         * OPTIONS is the SURROGATES Toolbox option structure given by
%           srgtsOptionSet, and SURROGATE is the respective surrogate
%           structure.
%
%     PRESSRMS = srgtsCrossValidation(OPTIONS, KFOLDS, IDXFOLDS):
%     does the computation of the cross-validation errors using k-fold
%     strategy. KFOLDS is the number of folds and IDXFOLDS is the vector of
%     indexes organized for NBPOINTS/k clusters. IDXFOLDS is obtained using
%     srgtsGetKfolds.
%
%     [PRESSRMS, eXV] = srgtsCrossValidation(...):
%     also returns a NBPOINTS-by-NBSURROGATES  matrix with the
%     cross-validation errors. Here the cross-validation errors obtained
%     when one data point is ignored and the surrogate is fitted to the
%     other (NBPOINTS - 1) points, with the procedure repeated for each
%     data point. In a point, the cross-validation error is:
%
%     eXV = Yhat - Y;
%
%     where Yhat is the value evaluated with the surrogate and Y is the
%     actual value of the function.
%
%     [PRESSRMS, eXV, yhatXV] = srgtsCrossValidation(...):
%     also returns the NBPOINTS-by-NBSURROGATES matrix of cross-validation
%     prediction.
%
%     [PRESSRMS, eXV, yhatXV, predvarXV] = srgtsCrossValidation(...):
%     also returns the NBPOINTS-by-NBSURROGATES matrix of cross-validation
%     prediction variance (prediction variance is implemented only for KRG
%     and PRS models; for all others it is going to be NaN).
%
%     [PRESSRMS, eXV, yhatXV, predvarXV, srgtSRGTXV] = srgtsCrossValidation(...):
%     also returns all SURROGATES created during cross-validation.
%
%     [PRESSRMS, eXV, yhatXV, predvarXV, srgtSRGTXV, ...
%     srgtOPTXV] = srgtsCrossValidation(...): also returns all
%     OPTION structures created during cross-validation.
%
%Example:
%     % basic information about the problem
%     myFN = @cos; % this could be any user-defined function
%     designspace = [0;     % lower bound
%                    2*pi]; % upper bound
%
%     % create DOE
%     NbPoints = 5;
%     X = linspace(designspace(1), designspace(2), NbPoints)';
%     Y = feval(myFN, X);
%
%     % fit surrogate models
%     srgtOPT = srgtsPRSSetOptions(X, Y);
%
%     % calculate cross validation errors, and PRESSRMS
%     [PRESSRMS, eXV] = srgtsCrossValidation(srgtOPT);
%
%     PRESSRMS =
%
%     1.1465
%
%     eXV =
%
%     1.5700
%    -0.8006
%     0.6003
%    -0.8006
%     1.5700
%
% Results may change for different "srgtOPT" structure.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Felipe A. C. Viana
% felipeacviana@gmail.com
% http://sites.google.com/site/felipeacviana
%
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

surro_predvar = PARAM.surro_predvar;

% check inputs
NbPoints = size(SrgtOPT.T,1);
if nargin <= 2
    kfolds   = NbPoints;
    idxFolds = (1:kfolds).';
end

if ~exist('SrgtSRGT','var')
    SrgtSRGT = [];
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run
PbKp = SrgtOPT.P;
TbKp = SrgtOPT.T;

NbPointsPerFold = NbPoints / kfolds; % which is also the number of clusters
Surr = upper(surro_predvar);%{'KRG','PRS'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create data structures
srgtSRGTXV = cell(kfolds, 1);
srgtOPTXV  = cell(kfolds, 1);
if kfolds ~= NbPoints
    for c1 = 1 : kfolds
        IDX(:,c1) = idxFolds(((c1 - 1)*NbPointsPerFold + 1) : c1*NbPointsPerFold);
    end
else
    IDX = idxFolds.';
end
predvarXVc  = cell(kfolds, 1);
yhatXVc     = cell(kfolds, 1);
ypredXVc    = cell(kfolds, 1);
parfor c2 = 1 : size(IDX,2)
    param = PARAM;
    srgtOPT = SrgtOPT;
    srgtSRGT = SrgtSRGT;
    surr = Surr;
    Idx = IDX;
    idx = Idx(:,c2);
%     yhatXVtmp = yhatXV;
%     predvarXVtmp = predvarXV;
    Pbkp = PbKp; Tbkp = TbKp;
    Ptest = Pbkp(idx,:);

    Ptraining = Pbkp; Ptraining(idx,:) = [];
    Ttraining = Tbkp; Ttraining(idx)   = [];

    if isempty(srgtSRGT)
        % new surrogates
%         eval(sprintf('surogtOPT = srgts%sSetOptions(Ptraining,Ttraining);', srgtOPT.SRGT));
%         eval(sprintf('surogtSRGT = srgts%sFit(surogtOPT);', srgtOPT.SRGT));
        [surogtSRGT,surogtOPT] = fit_surrogate(Ptraining,Ttraining,lower(srgtOPT.SRGT));
    else
        surogtOPT = srgtOPT;
        surogtSRGT = srgtSRGT;
    end

    srgtSRGTXV{c2} = surogtSRGT;
    srgtOPTXV{c2}  = surogtOPT;
%     eval(sprintf('yhatAux = %s_predict(Ptest,surogtSRGT);', lower(srgtOPT.SRGT)));
    yhatAux       = predict_funcval(Ptest,surogtSRGT,lower(srgtOPT.SRGT));
    yhatXVc{c2} = yhatAux;
    if (~strcmpi(param.predvar_type,'global') && param.strat == 1) || param.strat == 2
        predvarXVtmp = NaN(1,numel(surr));
        ypredXVtmp = NaN(1,numel(surr));
        for s = 1 : numel(surr)
            switch surr{s}
                case 'GP'
                    GPmodel = fit_surrogate(Ptraining,Ttraining,lower(surr{s}));
                    [ypredXVtmp(:,s),predvarXVtmp(:,s)] = gp_predict(Ptest,GPmodel);
                case {'KRIG','KRG'}
                    KRGmodel = fit_surrogate(Ptraining,Ttraining,lower(surr{s}));
                    [ypredXVtmp(:,s),predvarXVtmp(:,s)] = krg_predict(Ptest,KRGmodel);
                case 'PRS'
                    PRSmodel = fit_surrogate(Ptraining,Ttraining,lower(surr{s}));
                    [ypredXVtmp(:,s),predvarXVtmp(:,s)] = prs_predict(Ptest,PRSmodel);
            end
        end
        predvarXVc{c2} = predvarXVtmp;
        ypredXVc{c2} = ypredXVtmp;
    end
end
% cross-validation errors and PRESSRMS
yhatXV = cell2mat(yhatXVc);
eXV      = yhatXV - TbKp;
PRESSRMS = sqrt(mean(eXV.^2));

Rvec = NaN(1,numel(Surr));
model_predvar = NaN;
option_predvar = NaN;
if (~strcmpi(PARAM.predvar_type,'global') && PARAM.strat == 1) || PARAM.strat == 2
    ypredXV = cell2mat(ypredXVc);
    predvarXV = cell2mat(predvarXVc);
    errorXV = sqrt((ypredXV - repmat(TbKp(:),1,size(ypredXV,2))).^2);
    for r = 1 : size(errorXV,2)
        Rvec(:,r) = corr(errorXV(:,r),predvarXV(:,r),'type','kendall');
    end
    [~,maxCorrID] = max(Rvec);
    surroName = upper(Surr{maxCorrID});
    switch surroName
        case 'GP'
            [model_predvar,option_predvar] = fit_surrogate(PbKp,TbKp,surroName);
        case {'KRIG','KRG'}
            [model_predvar,option_predvar] = fit_surrogate(PbKp,TbKp,surroName);
        case 'PRS'
            [model_predvar,option_predvar] = fit_surrogate(PbKp,TbKp,surroName);
    end
end
return

function [model,option] = fit_surrogate(xtr,ytr,surroName)
surroName = lower(surroName);
switch surroName
    case 'gp'
        [model,option] = gp_train(xtr,ytr);
    case {'krig','krg'}
        [model,option] = krig_train(xtr,ytr);
    case 'prs'
        [model,option] = prs_train(xtr,ytr);
    case 'rbf'
        [model,option] = rbf_train(xtr,ytr);
    case 'rbnn'
        [model,option] = rbnn_train(xtr,ytr);
    case 'shep'
        [model,option] = shep_train(xtr,ytr);
    case 'svr'
        [model,option] = svr_train(xtr,ytr);
end
return

function [fpred,predstd] = predict_funcval(xtest,model,surroName)
switch surroName
    case 'gp'
        [fpred,predstd] = gp_predict(xtest,model);
    case {'krig','krg'}
        [fpred,predstd] = krig_predict(xtest,model);
    case 'prs'
        [fpred,predstd] = prs_predict(xtest,model);
    case 'rbf'
        fpred           = rbf_predict(xtest,model);
        predstd         = [];
    case 'rbnn'
        fpred           = rbnn_predict(xtest,model);
        predstd         = [];
    case 'shep'
        fpred           = shep_predict(xtest,model);
        predstd         = [];
    case 'svr'
        fpred           = svr_predict(xtest,model);
        predstd         = [];
end
return