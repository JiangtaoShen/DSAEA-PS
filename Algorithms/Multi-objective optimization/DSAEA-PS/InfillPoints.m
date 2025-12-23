function PopNew = InfillPoints(PopObj_1,PopObj_2,PopObj_3,PopDec_1,PopDec_2,PopDec_3,Problem,A1Obj,v,A1,mu,cigema)

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
[FN,~] = NDSort_CSDR(CPopObj,inf);
F1PopObj = CPopObjOrig(FN==1,:);
F1PopDec = CPopDecOrig(FN==1,:);
scoreIBEA = IGDValue(PopObj_1,F1PopObj);
scoreRVEA = IGDValue(PopObj_2,F1PopObj);
scoreCSDR = IGDValue(PopObj_3,F1PopObj);
[~,index] = min([scoreIBEA,scoreRVEA,scoreCSDR]);
if index == 1
    PopDec = PopDec_1;
    PopObj = PopObj_1;
elseif index == 2
    PopDec = PopDec_2;
    PopObj = PopObj_2;
else
    PopDec = PopDec_3;
    PopObj = PopObj_3;
end
PopDec = [PopDec;F1PopDec];
PopObj = [PopObj;F1PopObj];
[~,index_U] = unique(PopDec,'rows');
PopDec = PopDec(index_U,:);
PopObj = PopObj(index_U,:);
NVe = NoActive(A1Obj,v)
if ~mod(Problem.maxFE,ceil(Problem.FE*cigema)) && NVe ~= 0
    PopNew = Se_Div(PopDec,PopObj,A1Obj,v);
else
    PopNew = Se_IBEA(PopDec,PopObj,A1,mu);
end

end

