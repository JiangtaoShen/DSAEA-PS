function [f,g] = griewank(X)
%% if no input is given, return dimensions, bounds and minimum
if nargin == 0
    prob.nx = 20;
    prob.nf = 1;
    prob.ng = 0;
    for i=1:prob.nx
        prob.range(i,:) = [-600,600];
    end
    f = prob;
else
    X(X < -600) = -inf;   X(X > 600) = inf;
    % keep all values in the search domain
    if sum(sum(isinf(X))) > 0
        error('Griewank:RangeError','Search range out of bound');
    end
    
%% output function value
    [samples,dims] = size(X);
    
    I = 1:dims;
    XX = (X.^2)./4000;
    summs = sum(XX,2);
    theta = X./sqrt(repmat(I,samples,1));
    C = cos(theta);
    prods = prod(C,2);
    f = summs - prods + ones(samples,1);
    g = [];
end
end
