function PopNew = Se_Div(PopDec,PopObj,A1Obj,V)

%------------------------------- Copyright --------------------------------
% Copyright (c) 2021 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

A1Obj = A1Obj - repmat(min(A1Obj),size(A1Obj,1),1);
[NVe,Va]  = NoActive(A1Obj,V);
Ve = 1:size(V,1);
Ve(Va) = [];
Distance = pdist2(PopObj,PopObj);
Distance(logical(eye(length(Distance)))) = inf;
sDis = sort(Distance,2);
Div  = sDis(:,1) + 0.01*sDis(:,2);
if NVe ~= 0
    NCluster  = NVe;
    Ve        = V(Ve,:);
    [IDX,~]   = kmeans(Ve,NCluster);
    PopObj = PopObj - repmat(min(PopObj,[],1),size(PopObj,1),1);
    
    cosine = 1 - pdist2(Ve,Ve,'cosine');
    cosine(logical(eye(length(cosine)))) = 0;
    gamma  = min(acos(cosine),[],2);
    
    Angle  = acos(1-pdist2(PopObj,Ve,'cosine'));
    [~,associate] = min(Angle,[],2);
    APD_S  = ones(size(PopObj,1),1);
    for i = unique(associate)'
        current1 = find(associate==i);
        if ~isempty(current1)
            % Calculate the APD value of each solution
            APD = Angle(current1,i)/gamma(i);
            % Select the one with the minimum APD value
            APD_S(current1,:) = APD;
        end
    end
    
    Cindex = IDX(associate); % Solution to cluster
    Next = zeros(NCluster,1);
    
    for i = unique(Cindex)'
        
        solution_Best = [];
        current = find(Cindex==i);
        t = unique(associate(current));
        for j = 1:size(t,1)
            currentS = find(associate==t(j));
            [~,id] = min(APD_S(currentS,:),[],1);
            solution_Best = [solution_Best;currentS(id)];
        end
        [~,best] = min(APD_S(solution_Best,:),[],1);
        Next(i)     = solution_Best(best);
    end
    index  = Next(Next~=0);
    PopNew = PopDec(index,:);
else
    [~,index_2] = max(Div);
    PopNew = PopDec(index_2,:);
end


end