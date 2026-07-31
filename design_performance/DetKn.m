function detK = DetKn( X, sequence, theta,p,kernel,TENSORISED)
% function detK = DetKn( X, sequence, theta,p,kernel,TENSORISED)
%__________________________________________________________________________
% For X (d*n) an n-point design in [0,1]^d, returns the determinants detK of the 
% kernel matrices for the sequence of designs Xn(:,(1:i1)),...,Xn(:,(1:im)), with 
% [i_1 ... i_m] = sequence = a vector of m consecutive integers in {1,...,n}
% If sequence = k, detK is just the determinant of the kernel matrix for Xn(:,(1:k))
% The kernel used is specified by theta,p,kernel,TENSORISED (see calcR_switch.m)
% Requires the kernel to be translation invariant with K(x,x)=1 for all x
% Uses calcR_switch.m, calcrx_switch.m
%-----
% Author: L. Pronzato, 2022 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE.

smax=sequence(end);
[d,nmax]=size(X);
if smax>nmax
    display('sequence is too long for the nb. of design points Xn') 
    detK=NaN(size(sequence));
    return
end

LdetK=length(sequence);
detK=NaN(1,LdetK);
timebar= waitbar(0,'det-Kn...');

n=sequence(1);
Xn=X(:,1:n);
Kn=calcR_switch(Xn,theta,p,kernel,TENSORISED);
Knm1=inv(Kn);
    Knm1_all=NaN(nmax,nmax);
    Knm1_all(1:n,1:n)=Knm1;
detK(1)=det(Kn);

for i=1:LdetK-1
    waitbar(i/LdetK,timebar);
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
    detK(i+1)=detK(i)*tn;
end
close(timebar)
end

