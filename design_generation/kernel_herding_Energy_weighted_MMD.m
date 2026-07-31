function [ X,Time ] = kernel_herding_Energy_weighted_MMD( Xcand, nmax, theta,p,kernel,TENSORISED,...
                                                                tau, X_evaluate_P )
% function [ X,Time ] = kernel_herding_Energy_weighted_MMD( Xcand, nmax, theta,p,kernel,TENSORISED,...
%                                                               tau, X_evaluate_P )
%__________________________________________________________________________
% = Algorithm 6.3: Greedy energy minimization (when tau = 0)
% = Algorithm 6.4: Greedy minimization of a distorted squared MMD 
%   when tau~=0
% by kernel herding, see [Karvonen, Pronzato and Zhigljavsky, 2026], Chap. 6
% The design space is [0,1]^d.
% --> uses calcrx_switch.m and, if tau~=0 and TENSORISED=1, potential_energy.m 
%__________________________________________________________________________
% The design is searched in the "set" Xcand (d*Q) of candidate points in [0,1]^d
% When tau~=0, the potential of a measure mu for the kernel K is required
% (K is specified by kernel, theta and p, see calcrx_switch.m) 
% When X_evaluate_P is empty, we assume that mu is uniform in [0,1]^d, K is
%   tensorised (i.e., TENSORISED=1), and the potential is calculated by
%   potential_energy.m; no computation of big interdistances matrix is then involved.
% When a matrix X_evaluate_P (d*Q2) is provided, the potential is evaluated by
% the mean of K(x,zi) for zi in X_evaluate_P (==> a matrix Q2*Q is computed)
% When tau=NaN, we let tau increase with k, from 0.5 to 1, with
%   tau=((k-1)/(nmax-1)+1)/2;
%----- input variables
%   Xcand (d*Q) = a set of Q>=nmax candidate points in [0,1]^d 
%   nmax = number of points generated 
%   theta = inverse of correlation length of K
%   p = other scalar parameter in K (e.g., power of a Riesz kernel)
%   kernel = e.g., 'Matern32', or 'Riesz', or 'Gaussian', see calcR_switch.m 
%   TENSORISED = 1 for a tensorised (tensor-product) kernel
%   tau = scalar, 0 for energy minimization, 1 for MMD minimization, other
%       values for distorted squared MMD minimization
%       (if tau = 0, the first point X(:,1) is the point in Xcand closest to the 
%       center of the cube [0,1]^d)
%   X_evaluate_P (d*Q2) = a set of Q2 evaluation points zi to compute K(x,zi)
%       for all x in Xcand and hence potentials P(x) for all those x
%       (only required if tau ~= 0 and TENSORISED~=1)
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

t1=tic;
[d,Q]=size(Xcand);
Time=NaN(1,nmax); 

if Q<nmax
    disp('not enough candidate points, or nmax too large')
    X=NaN(d,nmax);
    return
end 
tau0=tau;

if tau==0
    % this is an energy minimisation problem: we don't need the potential
    % of mu for K, Pot_mu
    Pot_mu=zeros(1,Q);
else
    % Compute the potential Pot_mu=E{K(Xi,X)} for all Xi in X --> row vector 1*Q
    % once for all
    if TENSORISED~=1
        if isempty(X_evaluate_P)
            % exact calculation (for tensorised kernel only!)
            disp('a design set X_evaluate_P must be provided')
            X=NaN(d,nmax);
            return
        else
            % evaluate Pot_mu by discrete summation over X_evaluate_P
            Pot_mu=mean(calcrx_switch(X_evaluate_P,theta,p,Xcand,kernel,TENSORISED),1);
            if any(isinf(Pot_mu))
                disp('Xcand and X_evaluate_P must be disjoint')
                X=NaN(d,nmax);
                return
            end
        end
    else
        [ Pot_mu, ~ ] = potential_energy( Xcand, theta,p,kernel,1);
        if any(isinf(Pot_mu))
            disp('Bad choice of kernel: energy must be finite to have finite potential a.e.')
            X=NaN(d,nmax);
            return
        end
        
    end
end

% Initialisation
Xall=NaN(d,nmax);
IndX=zeros(1,nmax);

% Initialization valid whatever tau
if tau==0 || isnan(tau)
    % initialize at "central point" (in case of a design [0,1]^d) 
    Center=ones(d,1)/2; 
    D_to_center=pdist2(Center',Xcand');
    [~,i1]=min(D_to_center); 
else
    [~,i1]=max(tau*Pot_mu);
end    
X=Xcand(:,i1); IndX(1)=i1;
Xall(:,1)=X;
k=1;

% Initialize Sk=(1/k) sum_j=1^k K(zi,x_j) for all zi in Xcand --> will be a row vector 1*Q
Sk=0;

timebar= waitbar(0,'Kernel herding...'); 
while k<nmax
    waitbar(k/nmax,timebar);
    % vertex direction, step-size 1/(k+1)
    
    % Update Sk
    LastXj=X(:,end);
    Sk_update=calcrx_switch(LastXj,theta,p,Xcand,kernel,TENSORISED);
    Sk=(k-1)/k*Sk+Sk_update/k;

    if isnan(tau0)
        tau=((k-1)/(nmax-1)+1)/2;
    end
    FF=Sk-tau*Pot_mu;
    % avoid repetitions
        Ind_free=setdiff(1:Q,IndX);
        Dum=FF;
        Dum=Dum(Ind_free);
        [~,ibest]=min(Dum); ibest=ibest(1);
        imin=Ind_free(ibest);
      
    xnew=Xcand(:,imin);
    k=k+1;
    Xall(:,k)=xnew;
    X=Xall(:,1:k);
    IndX(k)=imin;  
    Time(k-1)=toc(t1);
end
Time(nmax)=toc(t1);
close(timebar)

end