function [ X,Time ] = greedy_MES( X_cand, nmax, theta,p,kernel,TENSORISED)
% function [ X,Time ] = greedy_MES( X_cand, nmax, theta,p,kernel,TENSORISED)
%__________________________________________________________________________
% = Algorithm 6.7: Greedy Maximum Entropy Sampling
%   see [Karvonen, Pronzato and Zhigljavsky, 2026], Chap. 6
% The design space is [0,1]^d.
% --> uses calcrx_switch.m 
%__________________________________________________________________________
% The design is searched in the "set" Xcand (d*N) of candidate points in [0,1]^d
% The kernel K is specified by kernel, theta and p, see calcrx_switch.m 
% Only applies to translation invariant kernels, with K(x,x)=1 for all x
%----- input variables
%   Xcand (d*N) = a set of N>=nmax candidate points in [0,1]^d 
%   nmax = number of points generated 
%   theta = inverse of correlation length of K
%   p = other scalar parameter in K (e.g., power of a Riesz kernel)
%   kernel = e.g., 'Matern32', or 'Riesz', or 'Gaussian', see calcR_switch.m 
%   TENSORISED = 1 for a tensorised (tensor-product) kernel
% The first point X(:,1) is always the closest one in Xcand to ones(d,1)/2 (i.e., the 
%   closest to the center of the cube [0,1]^d 
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

% first point (n=1) close to the center
dist_to_center=pdist2(0.5*ones(1,d),X_cand');
[~,i_close_to_center]=min(dist_to_center);
X=X_cand(:,i_close_to_center);
    Xall=NaN(d,nmax);
    Xall(:,1)=X;

Knm1=1; % --> we shall update Knm1=inv(Kn)
    Knm1_all=NaN(nmax,nmax);
    Knm1_all(1,1)=Knm1;
kn_C=calcrx_switch(X,theta,p,X_cand,kernel,TENSORISED); % now a vector 1*C, 
    kn_C_all=NaN(nmax,C);
    kn_C_all(1,:)=kn_C;
MES_crit=1-kn_C.*kn_C; MES_crit=max(MES_crit,0);
Time(1)=toc(t1);
timebar= waitbar(0,'greedy MES...'); 
for n=1:nmax-1
    waitbar(n/nmax,timebar);
    % choose next point

    [~,inew]=max(MES_crit);
    Xnew=X_cand(:,inew);
    an=kn_C(:,inew);
        Xall(:,n+1)=Xnew;
        X=Xall(:,1:n+1);    
    
    % update Knm1
    Knm1an=Knm1*an;    
    tn=1-an'*Knm1an;
        Knm1_all(1:n+1,1:n+1)=[Knm1 + Knm1an*Knm1an'/tn -Knm1an/tn
                               -Knm1an'/tn 1/tn];
        Knm1=Knm1_all(1:n+1,1:n+1);    
    % update kn_C
    bnp1=calcrx_switch(Xnew,theta,p,X_cand,kernel,TENSORISED);
     
    MES_crit=MES_crit-(Knm1an'*kn_C-bnp1).^2/tn; MES_crit=max(MES_crit,0);

        kn_C_all(1:n+1,:)=[kn_C
                           bnp1];
        kn_C=kn_C_all(1:n+1,:);      
    Time(n+1)=toc(t1);
end      
close(timebar)

end

