function [integrant_cdf] = integrant_cdf_symmetric_beta(rho,t,r,alpha,delta,d)
% function [integrant_cdf] = integrant_cdf_symmetric_beta(rho,t,r,alpha,delta,d)
%__________________________________________________________________________
% computes the integrand for F = Pr{||X-y|| <= t} when ||y||=r and 
% X has the spherically symmetric beta distribution 
% with parameters defined from alpha,delta,d through parameters_beta.m
% rho is a vector 1*n
% See Section 7.2 of the book [Karvonen, Pronzato and Zhigljavsky, 
% Springer, 2026] 
% 
% Uses density_norm_X.m
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

integrant_cdf=zeros(size(rho));

s=(t^2-(-r+rho).^2)./(4*r*rho); % first parameter for beta inverse 
Irh0_1=find(s>=1); % here the incomplete beta gives one
integrant_cdf(Irh0_1)=density_norm_X(rho(Irh0_1),alpha,delta,d); 

Irh0_good=find(s>0 & s<1);
integrant_cdf(Irh0_good)=density_norm_X(rho(Irh0_good),alpha,delta,d).*betainc(s(Irh0_good),(d-1)/2,(d-1)/2);

end    
