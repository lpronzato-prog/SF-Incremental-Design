function [ X,Time ] = greedy_IMSPE( X_cand, nmax, theta,p,kernel,TENSORISED,...
                                                                X_evaluate_P)
% function [ X,Time ] = greedy_IMSPE( X_cand, nmax, theta,p,kernel,PLOT,TENSORISED,...
%                                                                X_evaluate_P)
%__________________________________________________________________________
% = Algorithm 6.8: Greedy minimization of the IMSPE
%   see [Karvonen, Pronzato and Zhigljavsky, 2026], Chap. 6
% The design space is [0,1]^d.
% --> uses calcrx_switch.m and, if X_evaluate_P==[], potential_K2.m and 
%       IMSPE_integrals.m
%__________________________________________________________________________
% The design is searched in the "set" Xcand (d*N) of candidate points in [0,1]^d
% The kernel K is specified by kernel, theta and p, see calcrx_switch.m 
% Only applies to translation invariant kernels, with K(x,x)=1 for all x
%----- input variables
%   Xcand (d*N) = a set of N >=nmax candidate points in [0,1]^d 
%   nmax = number of points generated 
%   theta = inverse of correlation length of K
%   p = other scalar parameter in K (e.g., power of a Riesz kernel)
%   kernel = e.g., 'Matern32', or 'Riesz', or 'Gaussian', see calcR_switch.m 
%   TENSORISED = 1 for a tensorised (tensor-product) kernel
%   X_evaluate_P (d*Q) = a set of Q evaluation points zi to compute K(x,zi)
%       for all x in Xcand (an N*Q kernel matrix is computed)
%       When X_evaluate_P=[], the kernel is tensorised and must be either 
%       'Matern12', 'Matern32', or 'Gaussian', IMSPE is then computed 
%       for the uniform measure on [0,1]^d;
%       otherwise IMSPE is for the uniform measure on X_evaluate_P  
%       (and N*Q must not be too big !)
% The kernel K must be translation invariant with K(x,x)=1 for all x
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

compute_integrals=0;
if isempty(X_evaluate_P)
    compute_integrals=1;
end    

t1=tic;
Time=NaN(1,nmax);  

% compute potential for X_cand and K^2
if compute_integrals==1
    PK2_X_cand=potential_K2( X_cand, theta,p,kernel,1);
else
    % evaluate Pot_{K^2,mu} by discrete summation over X_evaluate_P
    PK2_X_cand=mean( (calcrx_switch(X_evaluate_P,theta,p,X_cand,kernel,TENSORISED)).^2,1 );
end    

% first point (n=1): maximum potential
[~,i1]=max(PK2_X_cand);
X=X_cand(:,i1);
IndX=i1;
    Xall=NaN(d,nmax);
    Xall(:,1)=X;
    IndXall=NaN(1,nmax);
    IndXall(1)=IndX;

Knm1=1; % --> we shall update Knm1=inv(Kn)
    Knm1_all=NaN(nmax,nmax);
    Knm1_all(1,1)=Knm1;
kn_C=calcrx_switch(X,theta,p,X_cand,kernel,TENSORISED); % now a vector 1*C, 
    kn_C_all=NaN(nmax,C);
    kn_C_all(1,:)=kn_C;
if compute_integrals==1
    sn_C=IMSPE_integrals( X_cand,X, theta,p,kernel,TENSORISED); % now a vector 1*C,
else
    sn_C= mean( (calcrx_switch(X_evaluate_P,theta,p,X,kernel,TENSORISED)*ones(1,C)).*...
                 calcrx_switch(X_evaluate_P,theta,p,X_cand,kernel,TENSORISED), 1);
end    
    sn_C_all=NaN(nmax,C);
    sn_C_all(1,:)=sn_C;
MES_crit=1-kn_C.*kn_C; MES_crit=max(MES_crit,0);
Sn=PK2_X_cand(i1);
    Sn_all=NaN(nmax,nmax);
    Sn_all(1,1)=Sn;
DIMSPE=PK2_X_cand+PK2_X_cand(i1)*kn_C.*kn_C-2*kn_C.*sn_C;
IMSPE_crit=DIMSPE./MES_crit;
IMSPE_crit(i1)=-Inf; % to force selection of NEW points
Time(1)=toc(t1);
timebar= waitbar(0,'Greedy IMSPE...'); 
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
    
    % update kn_C, sn_C, etc.
    bnp1=calcrx_switch(Xnew,theta,p,X_cand,kernel,TENSORISED);
    Deltanp1 = Knm1an'*kn_C-bnp1;
    MES_crit=MES_crit-Deltanp1.^2/tn; MES_crit=max(MES_crit,0);
    if compute_integrals==1
        snp1=IMSPE_integrals( X_cand,Xnew, theta,p,kernel,TENSORISED); % now a vector 1*C,
    else
        snp1= mean( (calcrx_switch(X_evaluate_P,theta,p,Xnew,kernel,TENSORISED)*ones(1,C)).*...
                     calcrx_switch(X_evaluate_P,theta,p,X_cand,kernel,TENSORISED), 1);
    end    
    vnp1=Knm1an'*sn_C-snp1;
    deltanp1=Sn*Knm1an-sn_C(:,inew);
    Knm1delta=Knm1*deltanp1;
    DIMSPE=DIMSPE+(2/tn)*Deltanp1.*(Knm1delta'*kn_C-vnp1)+...
        (1/tn)^2*(Deltanp1.^2)*(Knm1an'*deltanp1-Knm1an'*sn_C(:,inew)+PK2_X_cand(inew));
    DIMSPE=max(DIMSPE,0);
    
        kn_C_all(1:n+1,:)=[kn_C
                           bnp1];
        kn_C=kn_C_all(1:n+1,:);
        Knm1_all(1:n+1,1:n+1)=[Knm1 + Knm1an*Knm1an'/tn -Knm1an/tn
                                -Knm1an'/tn 1/tn];
        Knm1=Knm1_all(1:n+1,1:n+1);
        Sn_all(1:n+1,1:n+1)=[Sn sn_C(:,inew)
                            sn_C(:,inew)' PK2_X_cand(inew)];
        Sn=Sn_all(1:n+1,1:n+1);
        sn_C_all(1:n+1,:)=[sn_C
                  snp1];
        sn_C=sn_C_all(1:n+1,:);
    
    IMSPE_crit=DIMSPE./MES_crit;
    IMSPE_crit(IndX)=-Inf; % to force selection of NEW points
    Time(n+1)=toc(t1);
end      
close(timebar)

end

