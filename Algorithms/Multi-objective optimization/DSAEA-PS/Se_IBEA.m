function PopNew = Se_IBEA(PopDec,PopObj,A1,mu)
% Kriging selection in K-RVEA

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Jiangtao Shen
A1Obj = A1.objs;
zmin = min([A1.objs;PopObj]);
zmax = max([A1.objs;PopObj]);
A1Obj = A1Obj - repmat(zmin,size(A1Obj,1),1);
PopObj = PopObj - repmat(zmin,size(PopObj,1),1);
range  = zmax - zmin;
if 0.05*max(range) < min(range)
    A1Obj = A1Obj./repmat(range,size(A1Obj,1),1);
    PopObj = PopObj./repmat(range,size(PopObj,1),1);
end
[FN,~] = NDSort_CSDR(A1Obj,inf);
[FN2,~] = NDSort_CSDR(PopObj,inf);
A1Obj = A1Obj(FN==1,:);
PopObj = PopObj(FN2==1,:);
PopDec = PopDec(FN2==1,:);
CObj = [PopObj;A1Obj];
kappa = 0.05;
n_Pop = size(PopObj,1);
n_A1 = size(A1Obj,1);
Next_1 = 1 : n_Pop;
Next_2 = 1 : n_A1;
Next = 1:(n_Pop + n_A1);
[Fitness,I,C] = CalFitness(CObj);
while length(Next_1) > mu
    [~,x]   = min(Fitness(Next));
    Fitness = Fitness + exp(-I(Next(x),:)/C(Next(x))/kappa);
    if x <= n_Pop
        Next_1(x) = [];
        n_Pop = n_Pop - 1;
    else
        Next_2(x - n_Pop) = [];
        n_A1 = n_A1 - 1;
    end
    Next = 1:(n_Pop + n_A1);
end
PopNew = PopDec(Next_1,:);
end