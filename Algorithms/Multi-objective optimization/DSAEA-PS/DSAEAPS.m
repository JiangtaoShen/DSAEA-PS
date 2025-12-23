classdef DSAEAPS < ALGORITHM
    % <multi/many> <real> <expensive>
    % DSAEA-PS
    % This algorithm is written by Jiangtao Shen   
    methods
        function main(Algorithm,Problem)
            %% Parameter
            wmax = 20;
            mu = 5;
            %% Generate the reference points and population
            [V0,~] = UniformPoint(Problem.N + 10,Problem.M);
            v0 = eye(Problem.M);
            NI    = 11*Problem.D-1;
            P     = UniformPoint(NI,Problem.D,'Latin');
            A2    = SOLUTION(repmat(Problem.upper-Problem.lower,NI,1).*P+repmat(Problem.lower,NI,1));
            A1    = A2;
            THETA = 5.*ones(Problem.M,Problem.D);
            theta = 5.*ones(1,Problem.D);
            KModel = cell(1,Problem.M);
            RModel = cell(1,Problem.M);
            %% Optimization
            while Algorithm.NotTerminated(A2)
                A1Obj = A1.objs;
                A1Dec = A1.decs;
                V = V0.*repmat(max(A1Obj,[],1)-min(A1Obj,[],1),size(V0,1),1);
                v = v0.*repmat(max(A1Obj,[],1)-min(A1Obj,[],1),size(v0,1),1);
                %% Construct dominance relation model
                zmin = min(A1Obj);
                zmax = max(A1Obj);
                DA1Obj = A1Obj;
                DA1Obj = DA1Obj - repmat(zmin,size(DA1Obj,1),1);
                range  = zmax - zmin;
                if 0.05*max(range) < min(range)
                    DA1Obj = DA1Obj./repmat(range,size(DA1Obj,1),1);
                end
                [FrontNo,~] = NDSort_CSDR(DA1Obj,length(A1));
                Dmodel = rbf_build(A1Dec,FrontNo');
                %% Construct normal models
                for i = 1 : Problem.M
                    [mS, mY]   = dsmerge(A1Dec, A1Obj(:,i));
                    dmodel     = dacefit(mS,mY,'regpoly0','corrgauss',THETA(i,:),1e-5.*ones(1,Problem.D),100.*ones(1,Problem.D));                    
                    KModel{i}   = dmodel;                    
                    THETA(i,:) = dmodel.theta;
                end
                PopDec_1 = A1Dec;
                PopDec_2 = A1Dec;
                PopDec_3 = A1Dec;
                [PopObj_1,PopObj_2,PopObj_3,PopDec_1,PopDec_2,PopDec_3] = MOEAK(PopDec_1,PopDec_2,PopDec_3,Problem,wmax,KModel,V);
                CPopObjOrig = [PopObj_1;PopObj_2;PopObj_3];
                CPopDecOrig = [PopDec_1;PopDec_2;PopDec_3];
                CPopObj = [PopObj_1;PopObj_2;PopObj_3];
                zmin = min(CPopObj);
                zmax = max(CPopObj);
                range = zmax - zmin;
                CPopObj = CPopObj - repmat(zmin,size(CPopObj,1),1);
                if 0.05*max(range) < min(range)
                    CPopObj = CPopObj./repmat(range,size(CPopObj,1),1);
                end
                [FN,maxFN] = NDSort_CSDR(CPopObj,inf);
                %% Predict the dominance front number
                PopDec = CPopDecOrig;
                PopObj = CPopObjOrig;
                NP = size(PopDec,1);
                FNO = inf.*ones(1,NP);
                for i = 1:NP
                    fno = rbf_predict(Dmodel, A1Dec,PopDec(i,:));
                    FNO(i) = fno;
                end
                [aa,bb] = kmeans(FNO',maxFN);
                for i = 1:maxFN
                    [~,index] = min(bb);
                    bb(index) = inf;
                    idx = find(aa == index);
                    FNO(idx) = i;
                end
                ss = [FN + FNO;abs(FN - FNO)];
                [aa,~] = NDSort(ss',inf);
                indexF1 = find(aa == 1);
                if length(indexF1)> mu
                    F1Dec = PopDec(indexF1,:);
                    F1Obj = PopObj(indexF1,:);
                    PopNew = Se_IBEA(PopDec,PopObj,A1,mu);
                else
                    PopNew = PopDec(indexF1,:);
                end
                %% re-evaluate infilled points and update A1 and A2
                New     = SOLUTION(PopNew);
                A1 = UpdataArchive(A1,New);
                A2 = A1;
                [FN,~] = NDSort(A2.objs,inf);
                A2 = A2(FN == 1);
            end
        end
    end
end
