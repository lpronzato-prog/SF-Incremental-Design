function [ X, Time ] = greedy_quantization( Xcand, X_evaluate, nmax, s)
% function [ X, Time ] = greedy_quantization( Xcand, X_evaluate, nmax, s)
%__________________________________________________________________________
% = Algorithm 6.6: Greedy minimization of L_s-mean quantization error
%   see [Karvonen, Pronzato and Zhigljavsky, 2026], Chap. 6
% The design space is [0,1]^d.
%__________________________________________________________________________
% The design is searched in the "set" Xcand (d*N) of candidate points in [0,1]^d
% The criterion is evaluated for mu uniform on the set X_evaluate (d*Q) 
% --- which may or may not not intersect Xcand
% The criterion (to be minimized) is 
% (1/Q) sum_{k=1}^Q min_{i=1,...,n} ||x_i-z_k||^s with s>0
%----- input variables
%   Xcand (d*N) = a set of N>=nmax candidate points in [0,1]^d 
%   nmax = number of points generated 
%   X_evaluate_P (d*Q) = a set of Q evaluation points 
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

[d,N]=size(Xcand);
[~,Q]=size(X_evaluate);

t1=tic;
Time=NaN(1,nmax);  

%------------------
% 1) Initialization
%------------------
Dgrids=pdist2(Xcand',X_evaluate');   % Big N*Q interdistance matrix
D_to_Xn=Inf*ones(1,Q);          % 1*Q vector of min distances to design points
D_to_Xn_ps=D_to_Xn;
Dgrids_ps=Dgrids.^(s);     % Big N*Q interdistance matrix ^(s)

X=[]; % no design point yet
X_all=NaN(d,nmax);
IndX=zeros(1,nmax);

%--------------------
% 2) Start iterations
%--------------------

n=0;
iter=1;
timebar= waitbar(0,'greedy quantization...');
while n<nmax
    waitbar(n/nmax,timebar);
    iter=iter+1;
    Crit=sum(min(ones(N,1)*D_to_Xn_ps,Dgrids_ps),2)/Q;
    %[~,ibest]=min(Crit);
    % avoid repetitions
        Ind_free=setdiff(1:N,IndX);
        Dum=Crit;
        Dum=Dum(Ind_free);
        [~,ibest]=min(Dum); ibest=ibest(1);
        ibest=Ind_free(ibest);

    x_new=Xcand(:,ibest);
        % add design point to X and plot
      
        n=n+1;
            X_all(:,n)=x_new;
            X=X_all(:,1:n);
            IndX(n)=ibest;  
        % Update 1*Q vector of min (distances to design points)^s
        D_to_Xn_ps=min(D_to_Xn_ps,Dgrids_ps(ibest,:)); 
    Time(n)=toc(t1);
end
close(timebar)
end        



