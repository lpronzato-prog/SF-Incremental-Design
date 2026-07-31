function [ IMSPE ] = IntegratedMSPE( Xn, sequence, theta,p,cov_type)
% function [ IMSPE ] = IntegratedMSPE( Xn, sequence, theta,p,cov_type)
%__________________________________________________________________________
% For Xn (d*n) an n-point design, returns 
% IMSPE of size m = length(sequence), a row vector containing
% the integrated MSPEs of the sequence of designs Xn(:,(1:i1)),...,Xn(:,(1:im)), with 
% [i_1 ... i_m] = sequence = a vector of m consecutive integers in
% {1,...,n}, assuming the variance is 1
% If sequence = k, IMSPE is just the integrated MSPE of Xn(:,(1:k))
% The kernel used is specified by theta,p,cov_type (see calcR_switch.m)
% Requires the kernel to be translation invariant with K(x,x)=1 for all x
% Only applies to tensorised kernels of type Matern12, Matern32, Gaussian.
% Uses calcrx_switch.m, potential_K2, IMSPE_integrals
% The integration measure mu is uniform on [0,1]^d
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE.

TENSORISED=1;
nmax=sequence(end);
BK=NaN(1,nmax); % IMSPE=1-BK;
timebar= waitbar(0,'integrated MSPE...'); 
PXn=potential_K2( Xn, theta,p,cov_type,TENSORISED);

X=Xn(:,1);
Knm1=1; % --> we shall update Knm1=inv(Kn)

Sn=PXn(1);
    Sn_all=NaN(nmax,nmax);
    Sn_all(1,1)=Sn;
BK(1)=Sn;

for n=1:nmax-1
    waitbar(n/nmax,timebar);
    % choose next point
    Xnp1=Xn(:,n+1);
    knp1=calcrx_switch(X,theta,p,Xnp1,cov_type,TENSORISED);
    snp1=IMSPE_integrals( X,Xnp1, theta,p,cov_type,TENSORISED); snp1=snp1(:);
    Knm1knp1=Knm1*knp1;
    rhonp1=1-knp1'*Knm1knp1;
    BK(n+1)=BK(n)+(Knm1knp1'*Sn*Knm1knp1+PXn(n+1)-2*Knm1knp1'*snp1)/rhonp1;
    
    X=Xn(:,1:n+1);  
    % update Knm1, Sn
    Knm1=[Knm1 + Knm1knp1*Knm1knp1'/rhonp1 -Knm1knp1/rhonp1
          -Knm1knp1'/rhonp1 1/rhonp1];  
        Sn_all(1:n+1,1:n+1)=[Sn snp1
                            snp1' PXn(n+1)];
        Sn=Sn_all(1:n+1,1:n+1);
end 
IMSPE=-BK(sequence)+1;
close(timebar) 
end

