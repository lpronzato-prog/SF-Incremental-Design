function [ Vm, Em, S2 ] = KarhunenL_dim1( accuracy,m0,N,W,theta,p,cov_type )
% function [ Vm, Em, S2 ] = KarhunenL_dim1( accuracy,m0,W,N,theta,p,cov_type )
%__________________________________________________________________________  
% --> uses calcR_switch.m
%__________________________________________________________________________
%   Karhunen-Loève kernel decomposition for a measure with weights W(i) on [0,1]
% Uses N grid points 0,1/(N-1), 2/(N-1), ... , 1
% each one receiving weight Wi (with W = (W1,...,WN)) 
% Keep m eigenpairs, with m=min(m0,m1) and m1 such that eigenvalues smaller
% than accuracy*lambda1 are ignored
% theta,p,cov_type = kernel characteristics, see calcR_switch.m
% Em contains mp=m eigenvalues (ordered by decreasing values)
% Columns of Vm (N*mp) = mp eigenfunctions at the N grid points 
% S2 = total sum of eigenvalues
%-----
% L. Pronzato, 2025
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
%-----

% Kernel decomposition 
Xgrid=(0:1:N-1)/(N-1);
BigK=calcR_switch(Xgrid,theta,p,cov_type,0);
% weighting measure on the diagonal    
W=diag(W); % = probability measure

W12=diag(diag(sqrt(W))); Wm12=diag(diag(W.^(-1/2))); 
[Vm,Em]=eig(W12*BigK*W12); % Here sum(diag(Em)) < 1
%sum(diag(Em))
Vm=Wm12*Vm; % Propertly normalized : orthogonality & sum_j=1^N [Vm_i(x_j)]^2 Wj = 1 for all i=1,...,m

Em0=diag(Em); [Em,Im]=sort(Em0,'descend'); 
% Only keep m terms
Ikeep=find(Em/Em(1)>accuracy);
m1=length(Ikeep);
%[m0 m1]
m=min(m1,m0);
Em=Em(1:m); %Lambda=diag(Em); % Lambda = diagonal matrix of eigenvalues Em

Im=Im(1:m); Vm=Vm(:,Im); % Columns of Vm = m eigenfunctions at grid points

% Total sum of eigenvalues (up to N)
S2=sum(Em0);

end

