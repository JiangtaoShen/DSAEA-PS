function mmodel = mlp_train(X, Y)
% MLPMODEL - To construct MLP model
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
	[m1,nx] = size(X);
	[m2,ny] = size(Y);
	if m1 ~= m2
		error('X and Y must have the same number of rows')
	end

	ny = size(Y, 2);

	%range = minmax(X');
	%neuron = [5, 1];
	actfunc = {'tansig' 'tansig'};
	%mmodel = newff(range, neuron, actfunc, 'trainlm');

	mmodel = newff(X', Y', 6, actfunc);
	mmodel.trainParam.epochs = 5000;
	% mmodel.trainParam.goal = 1.e-2;
	mmodel.trainParam.showWindow = 0;

	% avoid splitting data into train/validate/test sets
	% mmodel.divideFcn = '';
	mmodel = train(mmodel, X', Y');
    
end
