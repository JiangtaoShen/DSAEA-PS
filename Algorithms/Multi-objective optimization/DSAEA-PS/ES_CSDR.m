function index = ES_CSDR(PopObj,N)
% The environmental selection of NSGA-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    %% Non-dominated sorting by CSDR
    zmin = min(PopObj);
    zmax = max(PopObj);
    PopObj = PopObj - repmat(zmin,size(PopObj,1),1);
    range  = zmax - zmin;
    if 0.05*max(range) < min(range)
        PopObj = PopObj./repmat(range,size(PopObj,1),1);
    end
    [FrontNo,MaxFNo] = NDSort_CSDR(PopObj,N);
    Next = FrontNo < MaxFNo;   
    %% Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(PopObj,FrontNo);
    %% Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    [~,index] = find(Next == true);
end