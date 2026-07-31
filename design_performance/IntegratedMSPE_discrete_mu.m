function [ IMSPE ] = IntegratedMSPE_discrete_mu( Xn, sequence, theta,p,kernel,TENSORISED, Xtest, k_fracfact )
% function [ IMSPE ] = IntegratedMSPE_discrete_mu( Xn, sequence, theta,p,kernel,TENSORISED, Xtest, k_fracfact )
%__________________________________________________________________________
% For Xn (d*n) an n-point design, returns 
% IMSPE of size m = length(sequence), a row vector containing
% the integrated MSPEs of the sequence of designs Xn(:,(1:i1)),...,Xn(:,(1:im)), with 
% [i_1 ... i_m] = sequence = a vector of m consecutive integers in
% {1,...,n}, assuming the variance is 1
% If sequence = k, IMSPE is just the integrated MSPE of Xn(:,(1:k))
% The kernel used is specified by theta,p,kernel,TENSORISED (see calcR_switch.m)
% Requires the kernel to be translation invariant with K(x,x)=1 for all x
% Uses calcrx_switch.m
% The integration measure mu is uniform on Xtest,
% a d*Q matrix, for instance of scrambled Sobol' points when the space is [0,1]^d:
%   pSmm = sobolset(d); pSmm = scramble(pSmm,'MatousekAffineOwen');
%   Xtest=(net(pS,2^19))';
% When k_fracfact=k>0, Xtest is completed by a k^d fractional factorial design 
%      (==> k should be kept small form large d, k=2 is reasonable)
% ---> one may use [ Xtest ] = candidate_set( d,q,k,center,Scramble )
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE.

[d,~]=size(Xn);

if isempty(Xtest)==0 && k_fracfact>0
    % complete by a k^d fractional factorial design
    k=k_fracfact;
    Complement=((fullfact(k*ones(1,d))-1)/(k-1))';
    Xtest=[Xtest Complement];
end
[~,N]=size(Xtest);

nmax=sequence(end);
BK=NaN(1,nmax); % IMSPE=1-BK;
timebar= waitbar(0,'integrated MSPE...'); 

X=Xn(:,1);
Kn=1;
Knm1=1; % --> we shall update Knm1=inv(Kn)
kn_N=calcrx_switch(X,theta,p,Xtest,kernel,TENSORISED);
    kn_N_all=NaN(nmax,N);
    kn_N_all(1,:)=kn_N; % now a vector 1*N

BK(1)=mean(kn_N.^2);

for n=1:nmax-1
    waitbar(n/nmax,timebar);
    % choose next point
    Xnp1=Xn(:,n+1);
    knp1=calcrx_switch(X,theta,p,Xnp1,kernel,TENSORISED);
    unp1=Knm1*knp1;
    rhonp1=1-knp1'*unp1;

    knp1_N=calcrx_switch(Xnp1,theta,p,Xtest,kernel,TENSORISED);
    BK(n+1)=BK(n)+ mean(sum((unp1'*kn_N- knp1_N ).^2,1))/rhonp1;

    X=Xn(:,1:n+1);  
    % update Knm1, kn_N_all
    Knm1=[Knm1 + unp1*unp1'/rhonp1 -unp1/rhonp1
          -unp1'/rhonp1 1/rhonp1];  

    kn_N_all(1:n+1,:)=[kn_N
                       knp1_N];
    kn_N=kn_N_all(1:n+1,:);      
end 
IMSPE=-BK(sequence)+1;
close(timebar) 
end

