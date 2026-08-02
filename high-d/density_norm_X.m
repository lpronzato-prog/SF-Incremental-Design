function [p_dR] = density_norm_X(rho,alpha,delta,d)
% function [p_dR] = density_norm_X(rho,alpha,delta,d)
%__________________________________________________________________________
% Computes the density of ||X|| at ||X||=rho when X has the spherically 
% symmetric beta distribution with parameters defined from 
% alpha,delta,d through parameters_beta.m
% rho is a vector 1*n
% See Section 7.1.6 of the book [Karvonen, Pronzato and Zhigljavsky, 
% Springer, 2026] 
% Uses parameters_beta.m
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

[a,b,M] = parameters_beta(alpha,delta,d);
% select rho in [0,sqrt(M)]
Irh0_good=find(rho>=0 & rho<=sqrt(M)); 
rho_good=rho(Irh0_good);
%
p_dR=zeros(size(rho));

% compute Beta
Beta=@(t) t.^(a-1).*(M-t).^(b-1)/(M^(a+b-1)*beta(a,b));
p_dR(Irh0_good) = 2*rho_good.*Beta(rho_good.^2);

end