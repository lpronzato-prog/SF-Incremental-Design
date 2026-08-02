function Edist = Edistortion_cube(alpha,delta,d,n,s)
% function Edist = Edistortion_cube(alpha,delta,d,n,s)
%__________________________________________________________________________
% Computes an approximation of the expected distortion ( = (Ls-mean quantization error)^s) 
% for a random design with n points independently distributed with the 
% spherically symmetric beta distribution computed by density_norm_X(rho,alpha,delta,d)
% See Section 7.2 of the book [Karvonen, Pronzato and Zhigljavsky, 
% Springer, 2026] 
%
% Uses parameters_beta.m, expected_cdf_symmetric_beta.m
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

[~,~,M1] = parameters_beta(alpha,delta,d);
[~,~,M2] = parameters_beta(1,1,d);
bound=sqrt(M1)+sqrt(M2);

% E{quantization error^s}
Edist = s*integral(@(r) r.^(s-1).*(-expected_cdf_symmetric_beta(r,n,alpha,delta,d)+1),0,bound);
end
