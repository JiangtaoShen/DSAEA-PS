function demo_krigpls
% clear;close all;
rng('default');
prob_name = 'griewank'; % test problems included here - ackley, rastrigin, rosenbrock
prob = feval(prob_name);

numpop = 400;
numtest = ceil(0.2*numpop);

% prob.nx = 15;
% prob.range(:,1) = zeros(prob.nx*2,1);
% prob.range(:,2) = ones(prob.nx*2,1);
% prob = feval(prob_name);
% prob.nx = 15;

% Nnoise = prob.nx;
% x = rand(numpop,prob.nx);
% x = pretreat(x,'center');
% [~,s,v] = svd(x);
% d = diag(s);
% d = d(1:10)/sum(d);
% T = x*v(:,1:10);
% X0 = T*v(:,1:10)';
% noise = randn(numpop,prob.nx+Nnoise);
% noise = 0.005*noise/max(max(abs(noise)));
% X = [X0 rand(numpop,Nnoise)]+noise; 
% y = T*(10:-1:1)';

X = bsxfun(@plus,prob.range(1:prob.nx,1)',bsxfun(@times,(prob.range(1:prob.nx,2)'),lhsdesign(numpop,prob.nx))); %-prob.range(1:prob.nx,1)'
% X = bsxfun(@plus,prob.range(:,1)',bsxfun(@times,(prob.range(:,2)'-prob.range(:,1)'),X));
y = feval(prob_name,X);
% y = y.*f;
% y = norm_dim(y,minmax(y'));
% y = norm_dim(f,minmax(f'));

if prob.nx == 2
    x = gridsamp(prob.range',50);
    y = feval(prob_name,x);
    x1 = reshape(x(:,1),50,50);
    x2 = reshape(x(:,2),50,50);
    f = reshape(y,50,50);
    surf(x1,x2,f);
end

[train_ids,test_ids] = find_medoids(X,numtest,prob.range);

xtrain = X(train_ids,:);
ytrain = y(train_ids,1);
xtest = X(test_ids,:);
ytest = y(test_ids,1);

disp('==========================================================');

disp('Model Building With PLS Started...');
tic
surro_dr = krigpls_train(xtrain,ytrain);  
t1 = toc;
ypred_w_ecr = krigpls_predict(xtest,surro_dr);
disp('Model Building With PLS Finished.');

disp('==========================================================');

disp('Model Building Without PLS (Full Model) Started...');
tic
surro_wo_dr = krig_train(xtrain,ytrain);  
t2 = toc;
ypred_wo_ecr = krig_predict(xtest,surro_wo_dr);
disp('Model Building Without PLS (Full Model) Finished.');

disp('==========================================================');

RMSE_w_dr = mean(sqrt((ytest-ypred_w_ecr).^2));
RMSE_wo_dr = mean(sqrt((ytest-ypred_wo_ecr).^2));

fprintf('With PLS:\n');
fprintf('\t Optimal Number of Principal Components: %d,\n',surro_dr.model.OptPCs);
fprintf('\t Model Building Time: %6.4f, \n',t1);
fprintf('\t RMSE: %6.4f. \n',RMSE_w_dr);

disp('==========================================================');

fprintf('Without PLS (Full Model):\n');
fprintf('\t Model Building Time: %6.4f, \n',t2);
fprintf('\t RMSE: %6.4f. \n',RMSE_wo_dr);

figure(1);
plot(ytest,ypred_w_ecr,'.',ytest,ytest,'r.'); grid on;
legend({'Predicted','True'});
title(['With PLS, Num. of PCs:',num2str(surro_dr.model.OptPCs),' and RMSE: ' num2str(RMSE_w_dr)]);
figure(2);
plot(ytest,ypred_wo_ecr,'.',ytest,ytest,'r.'); grid on;
title(['Without PLS (Full Model), RMSE: ' num2str(RMSE_wo_dr)]);
legend({'Predicted','True'});
return

function [train_ids,test_ids] = find_medoids(x,n,range)
x = norm_dim(x,range);
[~,c] = kmeans(x,n,'emptyaction','singleton');
test_ids = knnsearch(x,c);
train_ids = setdiff((1:size(x,1))',test_ids(:));
return

function x_norm = norm_dim(x,range)
N = size(x,1);
if iscell(range)
    B = nan(size(x,2),2);
    for i = 1:numel(range)
        B(i,:) = range{i};
    end
else
    B = range;
end
lb = B(:,1); LB = repmat(lb',N,1);
ub = B(:,2); UB = repmat(ub',N,1);

x_norm = (x-LB)./(UB-LB);
return

