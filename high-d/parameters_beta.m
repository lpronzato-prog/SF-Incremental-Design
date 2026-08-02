function [a,b,M] = parameters_beta(alpha,delta,d)
% function [a,b,M] = parameters_beta(alpha,delta,d)
%__________________________________________________________________________
% Computes the parameters of a spherically symmetric beta distribution matching 
% the first three moments of ||X||^2 when X in R^d ~ product of one-dimensional
% symmetric beta distributions on [-delta, delta] with shape parameter alpha.
% See Section 7.1.6 of the book [Karvonen, Pronzato and Zhigljavsky, 
% Springer, 2026] 
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

a= d* ( 4*(d+2)*alpha^2+4*(4*d-1)*alpha+3*(5*d-4) )/...
    ( 8*alpha^2*(d+2)+8*alpha*(4*d-1)+6*d );
M = delta^2* ( 4*(d+2)*alpha^2+4*(4*d-1)*alpha+3*d )/...
    ( 3*(2*alpha+1)^2 );
b=a*( M*(2*alpha+1)/(d*delta^2) -1 );
end