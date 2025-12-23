function [f,g]= ackley(x)
if nargin == 0
		prob.nx = 50;
		prob.nf = 1;
		prob.ng = 0;
		for i = 1:prob.nx
            prob.range(i,:) =[-32.768, 32.768];
		end
		f = prob;
else
    [f,g] = ackley_true(x);
end

end

function [f,g] = ackley_true(X)
for i = 1:size(X,1)
    x = X(i,:);
    f(i,1) = -20 * exp(-0.2 * sqrt( sum(x .* x)/ numel(x))) ...
        - exp(sum(cos(2*pi*x)) / numel(x)) ...
        + 20 + 2.7182818284590452353602874713526625;
end
g = [];
end

