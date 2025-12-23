function r = krigecr_matern_3_2(theta,d)
% Matern 3/2 correlation function
%=================================
%           h   d
%   r_i = prod prod (1 + sqrt(3)*theta_l*m_il) * exp(-sqrt(3)*theta_l*m_il) ;  where, i = 1,...,n, 
%          l=1 i=1                                                                    l = 1,...,h 
%                                                                                 and j = 1,...,d
%          
%
%   m_il = |Wstar_l (x_i_j - x'_i_j)|
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
term1 = (1 + sqrt(3)*Theta*abs(d)); 
term2 = exp(-sqrt(3)*Theta*abs(d));
r = prod(term1.*term2,2);
return
