function [ IMSPE ] = IntegratedMSPE_BLM( Xn, sequence, m,m1,n1,w1,theta,p,kernel)
% function [ IMSPE ] = IntegratedMSPE_BLM( Xn, sequence, m,m1,n1,w1,theta,p,kernel)
%__________________________________________________________________________
% For Xn (d*n) an n-point design, returns 
% IMSPE of size m = length(sequence), a row vector containing
% the integrated MSPEs of the sequence of designs Xn(:,(1:i1)),...,Xn(:,(1:im)), with 
% [i_1 ... i_m] = sequence = a vector of m consecutive integers in
% {1,...,n}, assuming the variance is 1
% If sequence = k, IMSPE is just the integrated MSPE of Xn(:,(1:k))
% The kernel used is specified by theta,p,cov_type (see calcR_switch.m)
% Requires the kernel to be translation invariant with K(x,x)=1 for all x
% The kernel is tensorised and IMSPE is calculated from a Bayesian Linear Model, 
% using Tensorised_KarhunenL.m, see this function for the meaning of m,m1,n1,w1,
% The variance is 1
% Uses Tensorised_KarhunenL.m
% The integration measure mu is uniform on [0,1]^d
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE.

nmax=sequence(end);
IMSPE=NaN(1,nmax); 
timebar= waitbar(0,'integrated MSPE with BLM...'); 

% Karhunen-Loève decomposition:
[ LAMBDA,PHI ]=Tensorised_KarhunenL( Xn,m,m1,n1,w1,theta,p,kernel );
% PHI (n*m) contains the values of the m eigenfunctions at the n design points
% LAMBDA (1*m) contains the m eigenvalues

s2m=(-(PHI.^2)*LAMBDA'+1)'; % 1*n vector of variances
DL=diag(LAMBDA);
MBm1=DL; % inverse of MB (Bayesian information matrix), initialized at diag(LAMBDA)

for n=0:nmax-1
    waitbar(n/nmax,timebar);
    % choose next point
    vnp1=MBm1*PHI(n+1,:)';
    MBm1=MBm1-vnp1*vnp1'/(s2m(n+1)+PHI(n+1,:)*vnp1);
    IMSPE(n+1)=trace(MBm1);
end 
IMSPE=IMSPE(sequence);
close(timebar) 
end

