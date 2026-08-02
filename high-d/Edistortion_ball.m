function Edist = Edistortion_ball(R,d,n,s)
% function Edist = Edistortion_ball(R,d,n,s)
%__________________________________________________________________________
% Computes the expected distortion ( = (Ls-mean quantization error)^s) for 
% a random design with n points independently uniformly distributed in the 
% ball B_d(0,R) when the reference measure mu is uniform in B_d(0,1).
% See Section 7.1.4 of the book [Karvonen, Pronzato and Zhigljavsky, 
% Springer, 2026] -- No approximation
%
% Uses integral_P_R_ball.m
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

psiball=@(r,d) d*r.^(d-1); % density of ||U|| for U uniform in B_d(0,1)

% integral_P_R_ball(tt,R,rr,d) computes 
%   integral_P_R=@(t,R,r,d) integral(@(R) betainc(sigma(t,R,r),(d-1)/2,(d-1)/2).*phiball_design(R,d,sig),0,R);
%   for all t=T(i) and r=rr(i), i=1,...,length(T)
% where
% sigma=@(t,rho,r) min(max((t.^2-(-r+rho).^2)./(4*r*rho),0),1); % = first parameter for beta inverse (truncated to [0,1])
% phiball_design=@(r,d,R) d*r.^(d-1)/R^d; % = density of ||X|| for X uniform in B_d(0,R)

% E{quantization error^s}
Edist = s*integral2(@(t,r) t.^(s-1).*(1-integral_P_R_ball(t,R,r,d)).^n.*psiball(r,d),0,Inf,0,1);

end