function y = mlp_predict(X, mmodel)
% MLPPREDICT - To predict using MLP model
%
% Call
%    y = mlp_predict(X, rmodel)
%
% Input
% X      : Data Points
% rmodel : MLP model obtained using mlp_model
%
% Output:
% y   : Predicted response
%

	% Check arguments
	if nargin ~= 2
		error('mlp_predict requires 2 input arguments')
	end

	y = sim(mmodel, X')';
end
