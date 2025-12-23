x = lhsdesign(200,10);
f = sin(sum(x.^2,2));
x1 = rand(1000,10);
f1 = sin(sum(x1.^2,2));
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,f,size(x,2));
csum = cumsum(PCTVAR(2,:)'/sum(PCTVAR(2,:)));
ncomp = (find(csum >= 0.9999, 1));
meanx = mean(x,1);
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] = plsregress(x,f,ncomp);
yfit = [ones(size(x,1),1) x]*BETA;
error = f - yfit;
model = dace_train(XS,error);
yfitp = [ones(size(x1,1),1) x1]*BETA;
x1c = x1 - repmat(meanx,size(x1,1),1);
x1d = x1c*stats.W(:,1:ncomp);
errorp = dace_predict(x1d,model);
fp = errorp+yfitp;
error1 = fp - f1;
model1 = dace_train(x,f);
fp1 = dace_predict(x1,model1);
error2 = fp1 - f1;
sqrt(mean(error1.^2))
sqrt(mean(error2.^2))
corr(fp,f1,'type','kendall')
corr(fp1,f1,'type','kendall')
