function Lct = Lebesgue_ct( X, sequence, Xtest, k_fracfact, theta,p,kernel,TENSORISED)
% function Lct = Lebesgue_ct( X, sequence, Xtest, k_fracfact, theta,p,kernel,TENSORISED)
%__________________________________________________________________________
% For X (d*n) an n-point design in [0,1]^d, returns a row vector of Lebesgue 
% constants computed on Xtest for the the sequence of designs Xn(:,(1:i1)),...,Xn(:,(1:im)), with 
% [i_1 ... i_m] = sequence = a vector of m consecutive integers in {1,...,n}
% If sequence = k, Lct is just the Lebesgue constant for Xn(:,(1:k))
% The kernel used is specified by theta,p,kernel,TENSORISED (see calcR_switch.m)
% Requires the kernel to be translation invariant with K(x,x)=1 for all x
% Uses calcR_switch.m, calcrx_switch.m
% Xtest is a d*Q matrix, for instance of scrambled Sobol' points when the space is [0,1]^d:
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

smax=sequence(end);
[d,nmax]=size(X);
if smax>nmax
    display('sequence is too long for the nb. of design points Xn') 
    Lct=NaN(size(sequence));
    return
end

if isempty(Xtest)==0 && k_fracfact>0
    % complete by a k^d fractional factorial design
    k=k_fracfact;
    Complement=((fullfact(k*ones(1,d))-1)/(k-1))';
    Xtest=[Xtest Complement];
end
[d,C]=size(Xtest);

LLct=length(sequence);
Lct=NaN(1,LLct);
timebar= waitbar(0,'Lebesgue ct...');

n=sequence(1);
Xn=X(:,1:n);
Kn=calcR_switch(Xn,theta,p,kernel,TENSORISED);
Knm1=inv(Kn);
    Knm1_all=NaN(nmax,nmax);
    Knm1_all(1:n,1:n)=Knm1;
kn_C=calcrx_switch(Xn,theta,p,Xtest,kernel,TENSORISED);
    kn_C_all=NaN(nmax,C);
    kn_C_all(1:n,:)=kn_C;
norm1_wn=sum(abs(Kn\kn_C),1);
Lct(1)=max(norm1_wn);

for i=1:LLct-1
    waitbar(i/LLct,timebar);
    n=sequence(i);
    Xnew=X(:,n+1);
    Xn=X(:,1:n);
    an=calcrx_switch(Xn,theta,p,Xnew,kernel,TENSORISED);   
    % update Knm1
    Knm1an=Knm1*an;    
    tn=1-an'*Knm1an;
        Knm1_all(1:n+1,1:n+1)=[Knm1 + Knm1an*Knm1an'/tn -Knm1an/tn
                               -Knm1an'/tn 1/tn];
        Knm1=Knm1_all(1:n+1,1:n+1);    
    % update kn_C
    bnp1=calcrx_switch(Xnew,theta,p,Xtest,kernel,TENSORISED);
    
        kn_C_all(1:n+1,:)=[kn_C
                           bnp1];
        kn_C=kn_C_all(1:n+1,:);      
    norm1_wn=sum(abs(Knm1*kn_C),1);
    Lct(i+1)=max(norm1_wn);
end
close(timebar)
end

