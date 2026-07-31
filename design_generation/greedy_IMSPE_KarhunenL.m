function [ X,Time ] = greedy_IMSPE_KarhunenL( X_cand, nmax, m,m1,n1,w1,theta,p,kernel)
% function [ X,Time ] = greedy_IMSPE_KarhunenL( X_cand, nmax, m,m1,n1,w1,theta,p,kernel)
%__________________________________________________________________________
% = Algorithm 6.9: Greedy minimization of the IMSPE with spectral
%   decomposition; see [Karvonen, Pronzato and Zhigljavsky, 2026], Chap. 6
% The design space is [0,1]^d.
% --> uses calcrx_switch.m and Tensorised_KarhunenL.m
%__________________________________________________________________________
% The design is searched in the "set" Xcand (d*N) of candidate points in [0,1]^d
% The kernel K is specified by kernel, theta and p, see calcrx_switch.m 
% Only applies to translation invariant kernels, with K(x,x)=1 for all x
% The kernel is tensorised and the IMSPE is calculated from a truncated KL
% decomposition, using Tensorised_KarhunenL.m, see this function for the
% meaning of m,m1,n1,w1,
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

TENSORISED=1;
[d,C]=size(X_cand);

t1=tic;
Time=NaN(1,nmax); 
% Karhunen-Loève decomposition:
[ LAMBDA,PHI ]=Tensorised_KarhunenL( X_cand,m,m1,n1,w1,theta,p,kernel );
% PHI (n*m) contains the values of the m eigenfunctions at the n design points
% LAMBDA (1*m) contains the m eigenvalues
LAMBDA2=LAMBDA.^2;
m=length(LAMBDA);

% first point (n=1): 
[~,i1]=max((PHI.^2)*LAMBDA2');  
X=X_cand(:,i1);
IndX=i1;
    Xall=NaN(d,nmax);
    Xall(:,1)=X;
    IndXall=NaN(1,nmax);
    IndXall(1)=IndX;

Kn=1;
Knm1=1; % --> we shall update Knm1=inv(Kn)
    Knm1_all=NaN(nmax,nmax);
    Knm1_all(1,1)=Knm1;
kn_C=calcrx_switch(X,theta,p,X_cand,kernel,TENSORISED); % now a vector 1*C, later a matrix
    kn_C_all=NaN(nmax,C);
    kn_C_all(1,:)=kn_C;
MES_crit=1-kn_C.*kn_C; MES_crit=max(MES_crit,0);
PHIn=PHI(i1,:); % now a vector 1*m, later a matrix
    PHIn_all=NaN(nmax,m);
    PHIn_all(1,:)=PHIn;
T=PHIn'*kn_C;   % an m*C matrix, with elements ell,j phi_ell(x1)*K(x1,X(:,j))

DIMSPE= LAMBDA2*(T-PHI').^2; % numerator of criterion
IMSPE_crit=DIMSPE./MES_crit;
IMSPE_crit(i1)=-Inf; % to force selection of NEW points
Time(1)=toc(t1);
timebar= waitbar(0,'Greedy IMSPE (Karhunen Loève)...'); 
for n=1:nmax-1
    waitbar(n/nmax,timebar);
    % choose next point
    [~,inew]=max(IMSPE_crit);
    Xnew=X_cand(:,inew);
    an=kn_C(:,inew);
        Xall(:,n+1)=Xnew;
        X=Xall(:,1:n+1);
        IndXall(n+1)=inew;
        IndX=IndXall(1:n+1);
   
    % update Knm1
    Knm1an=Knm1*an;
    tn=1-an'*Knm1an;
    
    % update kn_C, T, etc.
    bnp1=calcrx_switch(Xnew,theta,p,X_cand,kernel,TENSORISED);
    Deltanp1 = Knm1an'*kn_C-bnp1;
    MES_crit=MES_crit-Deltanp1.^2/tn; MES_crit=max(MES_crit,0);
    %
    T=T+((Knm1an'*PHIn-PHI(inew,:))'*Deltanp1)/tn;   
    DIMSPE= LAMBDA2*(T-PHI').^2; % numerator of criterion
    DIMSPE=max(DIMSPE,0);
    
        kn_C_all(1:n+1,:)=[kn_C
                           bnp1];
        kn_C=kn_C_all(1:n+1,:);    
        Knm1_all(1:n+1,1:n+1)=[Knm1 + Knm1an*Knm1an'/tn -Knm1an/tn
                               -Knm1an'/tn 1/tn];
        Knm1=Knm1_all(1:n+1,1:n+1);
        PHIn_all(1:n+1,:)=[PHIn
                           PHI(inew,:)];
        PHIn=PHIn_all(1:n+1,:);      
    IMSPE_crit=DIMSPE./MES_crit;
    IMSPE_crit(IndX)=-Inf; % to force selection of NEW points
    Time(n+1)=toc(t1);
end      
close(timebar)

end
