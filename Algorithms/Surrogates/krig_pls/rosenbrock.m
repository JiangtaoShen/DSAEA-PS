function [f,g] = rosenbrock(x)
if nargin == 0
		prob.nx = 50;
		prob.nf = 1;
		prob.ng = 0;
		for i = 1:prob.nx
            prob.range(i,:) = [-5.0, 10.0];
		end
		f = prob;
else
    [f,g] = rosenbrock_true(x);
end    
end

function [f,g] = rosenbrock_true(X)
for i = 1:size(X,1)
    x = X(i,:);
    f(i,1) = sum(100*(x(1:end-1) .^2 - x(2:end)) .^ 2 + (x(1:end-1)-1) .^2);
end
g = [];
end