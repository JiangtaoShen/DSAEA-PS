function r = krigecr_corrgauss(theta,d)
% Gaussian correlation function
%===============================
%           h   d
%   r_i = prod prod exp(-theta_l * m_il^2) ;  where, i = 1,...,n, l = 1,...,h and j = 1,...,d
%          l=1 j=1
%
%   m_il = |Wstar_l (x_i_j ? x'_i_j)|
%
% If length(theta) = 1, then the model is isotropic, i.e. all  theta_j = theta
%
% Input:
% -------------------------------------------------------------------------------------------------
% theta :  parameters in the correlation function
% d     :  m*n matrix with differences between given data points weighted
%          with Wstar obtained from ECR for each Principal Components 
% Output:
% -------------------------------------------------------------------------------------------------
% r     :  correlation value for a particular Principal Component

Theta = theta(:).';
r = prod(exp(d.^2 * -Theta),2);
return
