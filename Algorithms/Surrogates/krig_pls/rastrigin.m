function [f,g,x] = rastrigin(x)
if nargin == 0
		prob.nx = 30;
		prob.nf = 1;
		prob.ng = 0;
		for i = 1:prob.nx
            prob.range(i,:) = [-20.0, 20.0];
		end
		f = prob;
else
    [f,g] = rastrigin_true(x);
end
return

function [f,g] = rastrigin_true(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% INPUT:
%
% xx = [x1, x2, ..., xd]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,m] = size(x);
% summ = 0;
% for i = 1:n
% 	xi = x(i,:);
% 	sum = sum + (xi.^2 - 10*cos(2*pi.*xi));
% end
f(:,1) = sum(((x.^2-10*cos(2*pi.*x))+10),2);
% f(:,2) = f(:,1); 
% f(:,2) = x(:,1).^2;
g = [];
return