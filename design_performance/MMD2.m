function [ mmd2 ] = MMD2( Xn, sequence, theta,p,kernel)
% function [ mmd2 ] = MMD2( Xn, sequence, theta,p,kernel)
%__________________________________________________________________________
% For Xn (d*n) an n-point design, returns 
% mmd2 of size m = length(sequence), a row vector containing
% the squared MMDs of the sequence of designs Xn(:,(1:i1)),...,Xn(:,(1:im)), with 
% [i_1 ... i_m] = sequence = a vector of m consecutive integers in {1,...,n}
% If sequence = k, mmd2 is just the squared MMD of Xn(:,(1:k))
% The kernel used is specified by theta,p,kernel (see calcR_switch.m)
% The MMD is computed for the tensorised form of the kernel
% Uses potential_energy.m, calcR_switch.m, calcrx_switch.m
%-----
% Author: L. Pronzato, 2022 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE.

nmax=sequence(end);
mmd2=NaN(1,nmax); 
% Compute potential for all Xn(:,i) and energy 
[All_Bn,energy]=potential_energy( Xn, theta,p,kernel,1);
% Initialise (n=1)
An=calcR_switch(Xn(:,1),theta,p,kernel,1);
Bn=All_Bn(1);
mmd2(1)=energy+An-2*Bn;
timebar= waitbar(0,'Computing MMD...'); 
for n=1:nmax-1
    waitbar(n/nmax,timebar);
    knp1=calcrx_switch(Xn(:,1:n),theta,p,Xn(:,n+1),kernel,1);
    Knp1np1=calcR_switch(Xn(:,n+1),theta,p,kernel,1);
    An=(n^2*An + 2*sum(knp1)+Knp1np1)/(n+1)^2;
    Bn=(n*Bn+All_Bn(n+1))/(n+1);
    mmd2(n+1)=energy+An-2*Bn;
end   
mmd2=mmd2(sequence);
close(timebar) 
end