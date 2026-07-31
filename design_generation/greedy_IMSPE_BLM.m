function [ X,Time ] = greedy_IMSPE_BLM( X_cand, nmax, m,m1,n1,w1,theta,p,kernel)
% function [ X,Time ] = greedy_IMSPE_BLM( X_cand, nmax, m,m1,n1,w1,theta,p,kernel)
%__________________________________________________________________________
% = Algorithm 6.10: Greedy minimization of the IMSPE based on a Bayesian 
%   Linear Model; see [Karvonen, Pronzato and Zhigljavsky, 2026], Chap. 6
% The design space is [0,1]^d.
% --> uses Tensorised_KarhunenL.m
%__________________________________________________________________________
% The design is searched in the "set" Xcand (d*N) of candidate points in [0,1]^d
% The kernel K is specified by kernel, theta and p, see calcrx_switch.m 
% Only applies to translation invariant kernels, with K(x,x)=1 for all x
% The kernel is tensorised and IMSPE is calculated from a Bayesian Linear Model, 
% using Tensorised_KarhunenL.m, see this function for the meaning of m,m1,n1,w1,
% The variance is 1
% The integration measure mu is uniform on [0,1]^d
%----- input variables
%   Xcand (d*N) = a set of N>=nmax candidate points in [0,1]^d 
%   nmax = number of points generated 
%   m,m1,n1,w1: see Tensorised_KarhunenL.m
%   theta = inverse of correlation length of K
%   p = other scalar parameter in K (e.g., power of a Riesz kernel)
%   kernel = e.g., 'Matern32', or 'Riesz', or 'Gaussian', see calcR_switch.m 
%----- output variables
% X (d*nmax) = the design generated
% Time (1*nmax) = running time along iterations    
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

[d,C]=size(X_cand);

t1=tic;
Time=NaN(1,nmax); 
% Karhunen-Loève decomposition:
[ LAMBDA,PHI ]=Tensorised_KarhunenL( X_cand,m,m1,n1,w1,theta,p,kernel );
% PHI (C*m) contains the values of the m eigenfunctions at the C candidate points
% LAMBDA (1*m) contains the m eigenvalues

s2m=(-(PHI.^2)*LAMBDA'+1)'; % 1*C vector of variances
DL=diag(LAMBDA);
Hm=PHI';    % m*C
Vn=DL*Hm; % an m*C matrix
Pn=sum(Vn.*Hm,1); % a 1*C vector
Qn=sum(Vn.*Vn,1); % a 1*C vector
IMSPE_crit=Qn./(s2m+Pn);
% X=[];
% IndX=[];
    Xall=NaN(d,nmax);
    IndXall=NaN(1,nmax);


timebar= waitbar(0,'Greedy IMSPE (Bayesian Linear Model)...'); 
for n=0:nmax-1
    waitbar(n/nmax,timebar);
    % choose next point    
    [~,inew]=max(IMSPE_crit);
    Xnew=X_cand(:,inew);
        Xall(:,n+1)=Xnew;
        X=Xall(:,1:n+1);
        IndXall(n+1)=inew;
        IndX=IndXall(1:n+1);    

    % update Vn, Pn, Qn
    dnp1=s2m(inew)+Pn(inew);
    vnp1Hm=Vn(:,inew)'*Hm;
    Pn=Pn-vnp1Hm.^2/dnp1;
    vnp1Vn=Vn(:,inew)'*Vn;
    Qn=Qn-2*(vnp1Hm.*vnp1Vn)/dnp1+vnp1Hm.^2*(Vn(:,inew)'*Vn(:,inew))/dnp1^2;
    Vn=Vn-Vn(:,inew)*vnp1Hm/dnp1;
    IMSPE_crit=Qn./(s2m+Pn);
    IMSPE_crit(IndX)=-Inf; % to force selection of NEW points
    Time(n+1)=toc(t1);
end
close(timebar)

end
