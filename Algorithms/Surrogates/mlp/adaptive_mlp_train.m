function mmodel = adaptive_mlp_train(X, Y)
% ADAPTIVE_MLP_MODEL - To construct Adpative MLP model
%
% Call
%    rmodel = mlp_model(X, Y)
%
% Input
% X  : Data Points X(i,:), i=1,...,m
% 
% Output
% rmodel : MLP Model
%

% Check arguments
	if nargin ~= 2
		error('mlp_model requires 2 input arguments')
	end

	% Check design points
	[m1, nx] = size(X);
	[m2, ny] = size(Y);
	assert(m1 == m2, 'X and Y must have the same number of rows.');

	ni = round(nx/2) : max(nx,4);

	actfunc = {'tansig' 'purelin'};
	tmp_model = cell(1, length(ni));
	perf = zeros(1, length(ni));
	for i = 1:length(ni);
		neuron = [i];
		tmp_model{i} = newff(X', Y', neuron, actfunc);
		tmp_model{i}.trainParam.epochs = 5000;
		tmp_model{i}.trainParam.showWindow = false;
		[tmp_model{i}, tr] = train(tmp_model{i}, X', Y');
		perf(i) = tr.perf(tr.epoch == tr.best_epoch);
	end

	[tmp, I] = min(perf);
	mmodel = tmp_model{I};
end
