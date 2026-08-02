function PT = integral_P_R_ball(tt,R,rr,d)
% function PT = integral_P_R_ball(tt,R,rr,d)
%__________________________________________________________________________
% computes int I_v((d-1)/2,(d-1)/2)*phi(rho) drho
% where rho ~ density d*r.^(d-1)/R^d and v=min(max((T.^2-(-r+rho).^2)./(4*r*rho),0),1);
% See Section 7.1.4 of the book [Karvonen, Pronzato and Zhigljavsky, 
% Springer, 2026]
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

sigma=@(t,rho,r) min(max((t.^2-(-r+rho).^2)./(4*r*rho),0),1); % first parameter for beta inverse (truncated to [0,1])
phiball_design=@(r,d,R) d*r.^(d-1)/R^d; % = density of ||X|| for X uniform in B_d(0,R)

integral_P_R=@(t,R,r,d) integral(@(rho) betainc(sigma(t,rho,r),(d-1)/2,(d-1)/2).*phiball_design(rho,d,R),0,R);

PT=NaN*tt;
LT=length(tt); 
for i=1:LT
    t=tt(i);
    r=rr(i);
    PT(i)=integral_P_R(t,R,r,d);
end    
end