function [F] = cdf_symmetric_beta(t,r,alpha,delta,d)
% function [F] = cdf_symmetric_beta(t,r,alpha,delta,d)
%__________________________________________________________________________
% Computes F = Pr{||X-y|| <= t} when ||y||=r and 
% X has the spherically symmetric beta distribution 
% with parameters defined from alpha,delta,d through parameters_beta.m
% r is a vector 1*n
% See Section 7.2 of the book [Karvonen, Pronzato and Zhigljavsky, 
% Springer, 2026] 
% 
% Uses density_norm_X.m, parameters_beta.m, integrant_cdf_symmetric_beta.m
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

n=length(r);
F=NaN(size(r));
f=@(x) density_norm_X(x,alpha,delta,d); % for the case r=||y||=0

[~,~,M] = parameters_beta(alpha,delta,d);
for i=1:n
    rr=r(i);
    if rr==0
        % simple case!
        F(i)=integral(@(x)f(x),0,t);
    else
        % general case
        if alpha==0
            % This is a special case too
            rho=delta*sqrt(d);
            s=(t^2-(-rr+rho)^2)/(4*rr*rho); % first parameter for beta inverse 
            if s>=1
                F(i)=1;
            elseif s<=0
                F(i)=0;
            else
                F(i)=betainc(s,(d-1)/2,(d-1)/2);
            end 
        else
            ff=@(x) integrant_cdf_symmetric_beta(x,t,rr,alpha,delta,d); % for the general case 
            F(i)=integral(@(x)ff(x),0,sqrt(M));
        end
    end
end    
end    
