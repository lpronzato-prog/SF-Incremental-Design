function [ X,Time ] = lazy_greedy_relaxedCRKmu( Xcand, nmax, theta,p,kernel,TENSORISED,X_evaluate )
% function [ X,Time ] = lazy_greedy_relaxedCRKmu( Xcand, nmax, theta,p,kernel,TENSORISED,X_evaluate )
%__________________________________________________________________________
% = Algorithm 6.5: Greedy minimization of the double relaxation of covering
%   see [Karvonen, Pronzato and Zhigljavsky, 2026], Chap. 6
% = lazy-greedy version of 
%   greedy_relaxed_CRKmu( Xcand, nmax, theta,p,kernel,TENSORISED,X_evaluate )
%   (using supermodularity of the criterion, the global variable ALPHA_g permits
%   to evaluate the acceleration compared to greedy_relaxed_CRKmu.m)
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

global ALPHA_g

[d,NC]=size(Xcand);
[~,Q]=size(X_evaluate);

t1=tic;
Time=NaN(1,nmax);  

%------------------
% 1) Initialization
%------------------
KNQ=calcrx_switch(Xcand,theta,p,X_evaluate,kernel,TENSORISED); % Big N*Q kernel matrix

X=[]; % no design point yet
X_all=NaN(d,nmax);
indices=zeros(1,nmax);

nPot_xin=zeros(1,Q); % corresponds to n*Pot_(K,xi_n)=sum_{i=1}^n K(x_i,z_k) for the current n and all z_k 

%--------------------
% 2) Start iterations
%--------------------

%for n=1:nmax
n=0;
iter=1;
timebar= waitbar(0,'lazy greedy relaxed covering...');
ALPHA=zeros(1,nmax); % to check the gain of lazy-greedy versus greedy
while n<nmax
    waitbar(n/nmax,timebar);
    iter=iter+1;
    if n<=1
        % not lazy yet: we need 2 iterations to compute bounds on increments
        Candidate_nPot_xin=ones(NC,1)*nPot_xin+KNQ;
        Crit=sum(Candidate_nPot_xin.^(-1),2)/Q;
        [minCrit,ibest]=min(Crit);

        if n==0; Crit1=minCrit; end
        if n==1; bounds_on_increments=Crit1-Crit; end  
    else
        phikbest=minCrit;
        % let's be lazy now
        n_candidates=NC; % NC candidates at the beginning
        alpha_k=0; % to check the gain of lazy-greedy versus greedy       
        while n_candidates>0
            alpha_k=alpha_k+1;
            % find best potential candidate 
            [~,itry]=max(bounds_on_increments);  
            % Compute true value of criterion for X(:,itry)
            Candidate_nPot_xin=nPot_xin+KNQ(itry,:);
            Crit_itry=sum(Candidate_nPot_xin.^(-1))/Q;
            phikbest=min(phikbest,Crit_itry);
            % update bounds_on_increments, etc
            Crit(itry)=Crit_itry;                          
            bounds_on_increments(itry)=minCrit-Crit_itry;
            if n_candidates==1
                n_candidates=0; % we had n_candidates=1 at previous pass, we can stop
            else
                n_candidates=length(find(minCrit-bounds_on_increments<phikbest));
            end
        end
        ALPHA(n)=alpha_k/NC;
    end

    n=n+1;        
    Admissible=setdiff(1:NC,indices(1:n-1));
    [minCrit,ibest]=min(Crit(Admissible));
    ibest=Admissible(ibest);
    
    x_new=Xcand(:,ibest);
    % add design point to X     
        indices(n)=ibest;
        
        X_all(:,n)=x_new;
        X=X_all(:,1:n);
        nPot_xin=nPot_xin+KNQ(ibest,:);
    Time(n)=toc(t1);
    
end
close(timebar)
ALPHA_g=ALPHA;
end        



