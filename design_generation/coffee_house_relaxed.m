function [ X,Time ] = coffee_house_relaxed( Xcand,nmax,Alpha,rand_init )
% function [ X,Time ] = coffee_house_relaxed( Xcand,nmax,Alpha,rand_init )
%__________________________________________________________________________
% = Algorithm 6.2: Relaxed greedy packing 
%       of [Karvonen, Pronzato and Zhigljavsky, 2026]
% The design space is [0,1]^d.  
%__________________________________________________________________________
% The design X is searched within Xcand (d*Q), a set of Q of candidate points
% No computation of Q*Q interdistances matrix is involved
% alpha = Alpha(n) at iteration n defines the amount of relaxation
% Alpha = ones(1,nmax) gives the standard greedy-packing 
%   (= coffee-house algorithm)
% ---> at iteration k, with x* and x_i respectively one of the points (in
% Xcand) furthest away from the current design and one of the closest design points to x*,
% the next point is (1-alpha)*x_i+alpha*x*
% The CR and PR efficiencies of X are at least min(Alpha)/2 for all n <= nmax; 
% see Corollary 6.9 in [Karvonen, Pronzato and Zhigljavsky, 2026].
%----- input variables
%   Xcand (d*Q) = a set of Q >= nmax candidate points in [0,1]^d 
%   nmax = number of points generated 
%   Alpha (1*nmax) = relaxations coefficients that set the level of relaxation
%   rand_init = a scalar in [0,1]: the first point x1= X(:,1) is random in 
%       a cube with edgelength rand_init and centered at (1/2,...1/2) 
%----- output variables
% X (d*nmax) = the design generated
% Time (1*nmax) = running time along iterations    
%-----
% Author: L. Pronzato, 2022 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

[d,~]=size(Xcand);
t1=tic;
Time=NaN(1,nmax);  

%------------------
% 1) Initialization
%------------------
Xall=NaN(d,nmax);
% initialize by a random "central point" (in [0,1]^d) 
X=rand_init*(rand(d,1)-0.5)+0.5*ones(d,1);

Xall(:,1)=X;
D_to_X=pdist2(Xcand',X');  % distances to design points

%--------------------
% 2) Start iterations
%--------------------
timebar= waitbar(0,'Coffee-house relaxed...'); 
Time(1)=toc(t1);

for n=2:nmax
    waitbar(n/nmax,timebar);
    % Find one of the furthest points to current design
    [~,ifurthest]=max(D_to_X);
         
    % Find one of the closest design points
    D_to_furthest=pdist2(X',Xcand(:,ifurthest)');
    [~,idesign]=min(D_to_furthest);
    
    alpha_n=Alpha(n);
    Xnew=(1-alpha_n)*X(:,idesign)+alpha_n*Xcand(:,ifurthest);

    D_to_X=min(D_to_X,pdist2(Xcand',Xnew')); % update all D_to_X     
    x_last=Xnew; 
    Xall(:,n)=x_last;
    X=Xall(:,1:n);
    Time(n)=toc(t1);
end
close(timebar) 
end
        

