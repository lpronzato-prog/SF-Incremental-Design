function [ IMSPE ] = IntegratedMSPE_KarhunenL( Xn, sequence, m,m1,n1,w1,theta,p,kernel)
% function [ IMSPE ] = IntegratedMSPE_KarhunenL( Xn, sequence, m,m1,n1,w1,theta,p,kernel)
%__________________________________________________________________________
% For Xn (d*n) an n-point design, returns 
% IMSPE of size m = length(sequence), a row vector containing
% the integrated MSPEs of the sequence of designs Xn(:,(1:i1)),...,Xn(:,(1:im)), with 
% [i_1 ... i_m] = sequence = a vector of m consecutive integers in
% {1,...,n}, assuming the variance is 1
% If sequence = k, IMSPE is just the integrated MSPE of Xn(:,(1:k))
% The kernel used is specified by theta,p,cov_type (see calcR_switch.m)
% Requires the kernel to be translation invariant with K(x,x)=1 for all x
% The kernel is tensorised and the IMSPE is calculated from a truncated KL
% decomposition, using Tensorised_KarhunenL.m, see this function for the
% meaning of m,m1,n1,w1,
% The variance is 1
% Uses calcrx_switch.m, Tensorised_KarhunenL.m
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
timebar= waitbar(0,'integrated MSPE with KL decomposition...'); 

% Karhunen-Loève decomposition:
[ LAMBDA,PHI ]=Tensorised_KarhunenL( Xn,m,m1,n1,w1,theta,p,kernel );
% PHI (n*m) contains the values of the m eigenfunctions at the n design points
% LAMBDA (1*m) contains the m eigenvalues
LAMBDA2=LAMBDA.^2;

X=Xn(:,1);
Knm1=1; % --> we shall update Knm1=inv(Kn)
BK(1)=LAMBDA2*(PHI(1,:).^2)'; % = sum_i lambda_i^2 phi_i^2(x1)

for n=1:nmax-1
    waitbar(n/nmax,timebar);
    % choose next point
    Xnp1=Xn(:,n+1);
    knp1=calcrx_switch(X,theta,p,Xnp1,kernel,TENSORISED);
    Knm1knp1=Knm1*knp1;
    rhonp1=1-knp1'*Knm1knp1;
    BK(n+1)=BK(n)+LAMBDA2*((Knm1knp1'*PHI(1:n,:)-PHI(n+1,:)).^2)'/rhonp1;
    
    X=Xn(:,1:n+1);  
    % update Knm1
    Knm1=[Knm1 + Knm1knp1*Knm1knp1'/rhonp1 -Knm1knp1/rhonp1
          -Knm1knp1'/rhonp1 1/rhonp1];  
end 
IMSPE=-BK(sequence)+1;
close(timebar) 
end

