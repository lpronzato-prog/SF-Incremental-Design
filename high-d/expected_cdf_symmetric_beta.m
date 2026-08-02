function [ECDFn] = expected_cdf_symmetric_beta(T,n,alpha,delta,d)
% function [ECDFn] = expected_cdf_symmetric_beta(T,n,alpha,delta,d)
%__________________________________________________________________________
% For all t in T, computes F = E_U{ Pr_Xn {min_i ||U-X_i||<=t} } 
% = 1 - E_U{ [1-Pr{ ||U-X_i|| <= t }]^n}
% when the n points X_i are iid with the spherically symmetric beta distribution 
% computed by density_norm_X(rho,alpha,delta,d) (with rho = ||X_i||) 
% and U has the spherically symmetric beta distribution with density density_norm_X(r,1,1,d)
% (with r = ||U||) 
% n = number of design points
% t is a vector 1*LT, returns a vector ECDFn of the same size
% See Section 7.2 of the book [Karvonen, Pronzato and Zhigljavsky, 
% Springer, 2026] 
% 
% Uses density_norm_X.m, parameters_beta.m, cdf_symmetric_beta.m
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE.

[~,~,M] = parameters_beta(1,1,d);
LT=length(T);
ECDFn=zeros(1,LT);
for i=1:LT
    t=T(i);
    cdf_given_r=@(r) cdf_symmetric_beta(t,r,alpha,delta,d);
    ECDFn(i)=1-integral(@(r) (1-cdf_given_r(r)).^n.*density_norm_X(r,1,1,d) ,0,sqrt(M));
end

end