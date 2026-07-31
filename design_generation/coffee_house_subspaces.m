function [ X,Time ] = coffee_house_subspaces( Xcand, nmax, dims, Wdim, randini, randdim)
% function [ X,Time ] = coffee_house_subspaces( Xcand, nmax, dims, Wdim, randini, randdim)
%__________________________________________________________________________
% = Algorithm 6.1-PR(D,W): Greedy packing with projections 
%       of [Karvonen, Pronzato and Zhigljavsky, 2026], Section 6.4.2
% The design space is [0,1]^d. 
% --> uses tri2p01.m 
%__________________________________________________________________________
% The design is searched in the "set" Xcand (d*Q) of candidate points in [0,1]^d
% No computation of Q*Q interdistances matrix is involved
%----- input variables
%   Xcand (d*Q) = a set of Q >= nmax candidate points in [0,1]^d 
%   nmax = number of points generated 
%   dims = dimensions of subspaces of interest (e.g., dims=[1 2] to consider 
%           projections in one and two dimensional subspaces)
%   Wdim = weights (importance) allocated to each dimension -- same size as dims
%   randini = 1 ==> the first point X(:,1) is randomized within Xcand
%           otherwise, X(:,1) is the closest one in Xcand to ones(d,1)/2 (i.e., the 
%           closest to the center of the cube [0,1]^d) 
%   randdim = 1 ==> the choice of next candidate is randomized in case of
%           several equivalent solutions
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

[d,Q]=size(Xcand);
t1=tic;
Time=NaN(1,nmax); 
dims=dims(:); 
MMproj=tri2p01(d); MMproj=MMproj(:,1:d); % indicator matrix for all possible projections                                          
iselectproj=ismember(sum(MMproj,2),dims);
MMproj=MMproj(iselectproj==1,:); % m*d matrix with elements in {0,1} defining the m projections of interest                              
[m,~]=size(MMproj);

%------------------
% 1) Initialization
%------------------
Xall=NaN(d,nmax); % to avoid using a matrix of increasing size
if randini==1
    i1=randi(Q,1); % random point in Xcand
else
    % initialize by grid point closest to "central point" (in case of a design [0,1]^d)
    Center=ones(d,1)/2; 
    D_to_center=pdist2(Center',Xcand');
    [~,i1]=min(D_to_center); 
end    
X=Xcand(:,i1); 
Xall(:,1)=X;
% Initialize the Q*m matrix of distances to design points for all m
%       projections (including the full space)
D_to_X_proj=NaN(Q,m); 
Weights=NaN(1,m);
for i=1:m
    iproj=find(MMproj(i,:)==1);
    D_to_X_proj(:,i)=pdist2(Xcand(iproj,:)',X(iproj)');
    Weights(i)=Wdim(dims==sum(MMproj(i,:))); % weight for that projection
end    

%--------------------
% 2) Start iterations
%--------------------
timebar= waitbar(0,'Coffee-house in subspaces...'); 
Time(1)=toc(t1);
for n=2:nmax
    waitbar(n/nmax,timebar);
    
    [Distmax,~]=max(D_to_X_proj);  % maximum distance from design (for each projection)
    D_to_X_proj_normalized=D_to_X_proj.*(ones(Q,1)*Distmax.^(-1)); % normalized distances to design
    Weighted_D_to_X_proj_normalized=D_to_X_proj_normalized.*(ones(Q,1)*(Weights).^(-1)); % weighted normalized distances to design
    Weighted_D_to_X_proj_normalized_min=min(Weighted_D_to_X_proj_normalized,[],2);   % minimum over projections of normalized distances 
    [distmaxmin,~]=max(Weighted_D_to_X_proj_normalized_min); % --> we want the point that maximizes this minimum
    
    Imax=find(Weighted_D_to_X_proj_normalized_min==distmaxmin);
    if length(Imax)>=1
        ibest=Imax(1);
        if randdim==1
            ibest=Imax(randi(length(Imax),1)); % randomisation!
        end
    end
    Xnew=Xcand(:,ibest);
    % update D_to_X_proj = distances to design for each subspace of interest
    D_to_X_proj_new=NaN(Q,m); 
    for i=1:m
        iproj=find(MMproj(i,:)==1);
        D_to_X_proj_new(:,i)=pdist2(Xcand(iproj,:)',Xnew(iproj)');
    end 
    D_to_X_proj=min(D_to_X_proj,D_to_X_proj_new); 
    
    x_last=Xnew; 
    Xall(:,n)=x_last;
    X=Xall(:,1:n);
     
    Time(n)=toc(t1);
end
close(timebar) 

end
        