function [f,g,x] = TF(x)
if nargin == 0
		prob.nx = 30;
		prob.nf = 1;
		prob.ng = 0;
		for i = 1:prob.nx
            prob.range(i,:) = [-20.0, 20.0];
		end
		f = prob;
else
    [f,g] = tf_true(x);
end
return

function [f,g] = tf_true(x)
% [~,m] = size(x);
f(:,1) = sum(((x.^2-10*cos(2*pi.*x))+10),2);
g = [];
return