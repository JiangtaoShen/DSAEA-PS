function [PopObj_1,PopObj_2, PopObj_3, PopDec_1, PopDec_2, PopDec_3] = MOEAK(PopDec_1,PopDec_2,PopDec_3,Problem,wmax,KModel,V)
w1      = 1;
w2      = 1;
w3      = 1;
% IBEA
while w1 <= wmax
    OffDec_1 = OperatorGA(PopDec_1);PopDec_1 = [PopDec_1;OffDec_1];
    [N,~]  = size(PopDec_1);
    PopObj_1 = zeros(N,Problem.M);
    for i = 1: N
        for j = 1 : Problem.M
            PopObj_1(i,j) = predictor(PopDec_1(i,:),KModel{j});
        end
    end
    index_1  = ES_IBEA(PopObj_1,Problem.N);PopDec_1 = PopDec_1(index_1,:);PopObj_1 = PopObj_1(index_1,:);
    w1 = w1 + 1;
end
% RVEA
while w2 <= wmax
    OffDec_2 = OperatorGA(PopDec_2);PopDec_2 = [PopDec_2;OffDec_2];
    [N,~]  = size(PopDec_2);
    PopObj_2 = zeros(N,Problem.M);
    for i = 1: N
        for j = 1 : Problem.M
            PopObj_2(i,j) = predictor(PopDec_2(i,:),KModel{j}); 
        end
    end
    index_2  = ES_RVEA(PopObj_2,V,(w2/wmax).^2);PopDec_2 = PopDec_2(index_2,:);PopObj_2 = PopObj_2(index_2,:);
    w2 = w2 + 1;
end
% NSGA-II/CSDR
while w3 <= wmax
    drawnow();
    OffDec_3 = OperatorGA(PopDec_3);PopDec_3 = [PopDec_3;OffDec_3];
    [N,~]  = size(PopDec_3);
    PopObj_3 = zeros(N,Problem.M);
    for i = 1: N
        for j = 1 : Problem.M
            PopObj_3(i,j) = predictor(PopDec_3(i,:),KModel{j});
        end
    end
    index_3  = ES_CSDR(PopObj_3,Problem.N);PopDec_3 = PopDec_3(index_3,:);PopObj_3 = PopObj_3(index_3,:);
    w3 = w3 + 1;
end
end

