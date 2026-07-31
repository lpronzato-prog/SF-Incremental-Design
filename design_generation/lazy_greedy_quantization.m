function [ X,Time ] = lazy_greedy_quantization( Xgrid, Xcand, nmax, s)
% function [ X,Time ] = lazy_greedy_quantization( Xgrid, Xcand, nmax, s)
%__________________________________________________________________________
% = Algorithm 6.6: Greedy minimization of L_s-mean quantization error
%   see [Karvonen, Pronzato and Zhigljavsky, 2026], Chap. 6
% = lazy-greedy version of greedy_quantization( Xcand, Xgrid, nmax, s)
%   (using supermodularity of the criterion, the global variable ALPHA_g permits
%   to evaluate the acceleration compared to greedy_quantization.m)
% The design space is [0,1]^d.
%__________________________________________________________________________
% The design is searched in the "set" Xcand (d*N) of candidate points in [0,1]^d
% The criterion is evaluated for mu uniform on the set Xgrid (d*Q) 
% --- which may or may not not intersect Xcand
% The criterion (to be minimized) is 
% (1/Q) sum_{k=1}^Q min_{i=1,...,n} ||x_i-z_k||^s with s>0
%----- input variables
%   Xcand (d*N) = a set of N>=nmax candidate points in [0,1]^d 
%   nmax = number of points generated 
%   Xgrid (d*Q) = a set of Q evaluation points 
%   (an N*Q matrix of interdistances is computed) 
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
[~,Q]=size(Xgrid);

t1=tic;
Time=NaN(1,nmax);  
powerq=s;

%------------------
% 1) Initialization
%------------------
Dgrids=pdist2(Xcand',Xgrid');   % Big NC*Q interdistance matrix
D_to_Xn=Inf*ones(1,Q);          % 1*Q vector of min distances to design points
D_to_Xn_pq=D_to_Xn;             % 1*Q vector of (min distances to design points)^q
Dgrids_pq=Dgrids.^(powerq);     % Big NC*Q interdistance matrix ^(q)

X=[]; % no design point yet
X_all=NaN(d,nmax);
indices=zeros(1,nmax);
%--------------------
% 2) Start iterations
%--------------------

%for n=1:nmax
n=0;
iter=1;
timebar= waitbar(0,'lazy greedy quantization...');
ALPHA=zeros(1,nmax); % to check the gain of lazy-greedy versus greedy
while n<nmax
    waitbar(n/nmax,timebar);
    iter=iter+1;
    if n<=1
        % not lazy yet: we need 2 iterations to compute bounds on increments
        Crit=sum(min(ones(NC,1)*D_to_Xn_pq,Dgrids_pq),2)/Q;
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
            Crit_itry=sum(min(D_to_Xn_pq,Dgrids_pq(itry,:)))/Q;
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
    %----
    ibest=ibest(1);
    %----
    x_new=Xcand(:,ibest);
    % add design point to X and plot     
   
        indices(n)=ibest;
        X_all(:,n)=x_new;
        X=X_all(:,1:n);
        D_to_Xn_pq=min(D_to_Xn_pq,Dgrids_pq(ibest,:));
    Time(n)=toc(t1);   
end
close(timebar)
ALPHA_g=ALPHA;
end


