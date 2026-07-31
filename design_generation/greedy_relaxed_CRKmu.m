function [ X, Time ] = greedy_relaxed_CRKmu( Xcand, nmax, theta,p,kernel,TENSORISED,X_evaluate )
% function [ X, Time ] = greedy_relaxed_CRKmu( Xcand, nmax, theta,p,kernel,TENSORISED,X_evaluate )
%__________________________________________________________________________
% = Algorithm 6.5: Greedy minimization of the double relaxation of covering
%   see [Karvonen, Pronzato and Zhigljavsky, 2026], Chap. 6
% The design space is [0,1]^d.
% --> uses calcrx_switch.m 
%__________________________________________________________________________
% The design is searched in the "set" Xcand (d*N) of candidate points in [0,1]^d
% The kernel K is specified by kernel, theta and p, see calcrx_switch.m 
% The criterion is evaluated for mu uniform on the set X_evaluate (d*Q) 
% --- which may or may not not intersect Xcand
% The relaxed criterion (to be minimized) is 
% (1/Q) sum_{k=1}^Q ( sum_{i=1}^n K(x_i,z_k) )^{-1} 
%----- input variables
%   Xcand (d*N) = a set of N>=nmax candidate points in [0,1]^d 
%   nmax = number of points generated 
%   theta = inverse of correlation length of K
%   p = other scalar parameter in K (e.g., power of a Riesz kernel)
%   kernel = e.g., 'Matern32', or 'Riesz', or 'Gaussian', see calcR_switch.m 
%   TENSORISED = 1 for a tensorised (tensor-product) kernel
%   X_evaluate_P (d*Q) = a set of Q evaluation points zi to compute K(x,zi)
%       for all x in Xcand (an N*Q kernel matrix is computed)
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

[d,N]=size(Xcand);
[~,Q]=size(X_evaluate);

t1=tic;
Time=NaN(1,nmax);  

%------------------
% 1) Initialization
%------------------
KNQ=calcrx_switch(Xcand,theta,p,X_evaluate,kernel,TENSORISED); % Big N*Q kernel matrix
% 
X=[]; % no design point yet
X_all=NaN(d,nmax);
IndX=zeros(1,nmax);

nPot_xin=zeros(1,Q); % corresponds to n*Pot_(K,xi_n)=sum_{i=1}^n K(x_i,z_k) for the current n and all z_k 
%--------------------
% 2) Start iterations
%--------------------

%for n=1:nmax
n=0;
iter=1;
timebar= waitbar(0,'greedy_relaxed_covering...');
while n<nmax
    waitbar(n/nmax,timebar);
    iter=iter+1;
        Candidate_nPot_xin=ones(N,1)*nPot_xin+KNQ;
        Crit=sum(Candidate_nPot_xin.^(-1),2)/Q;
    
    % avoid repetitions
        Ind_free=setdiff(1:N,IndX);
        Dum=Crit;
        Dum=Dum(Ind_free);
        [~,ibest]=min(Dum); ibest=ibest(1);
        ibest=Ind_free(ibest);

    x_new=Xcand(:,ibest);
        % add design point to X
        
        n=n+1;
        X_all(:,n)=x_new;
        X=X_all(:,1:n);
        IndX(n)=ibest;  
        nPot_xin=nPot_xin+KNQ(ibest,:);
    Time(n)=toc(t1);
end
close(timebar)
end        



