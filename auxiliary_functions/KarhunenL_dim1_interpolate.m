function [ PHIx ] = KarhunenL_dim1_interpolate( Vm,Em,W,theta,p,cov_type,x  )
% function [ PHIx ] = KarhunenL_dim1_interpolate( Vm,Em,W,theta,p,cov_type,x )
%__________________________________________________________________________  
% --> uses calcrx_switch.m
%__________________________________________________________________________
%   Based on the Karhunen-Loève kernel decomposition for a measure with weights W(i) on [0,1]
% Uses the decomposition Vm, Em at N grid points 0,1/(N-1), 2/(N-1), ... , 1
% obtained by KarhunenL_dim1.m
% Em = m eigenvalues (decreasing order)
% W = N weights of mu at the grid points
% Columns of Vm (N*m) = m eigenfunctions at the N grid points 
% theta,p,cov_type = kernel characteristics, see calcR_switch.m
% x = a (1*NN) row vector of locations where we want to compute the eigenfunctions
% PHIx (NN*m) = the m eigenfunctions at x 
%-----
% L. Pronzato, 2025
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
%-----

[N,m]=size(Vm);
Xgrid=(0:1:N-1)/(N-1);
W=diag(W);

% Original kernel
KxXgrid=(calcrx_switch(Xgrid,theta,p,x,cov_type,0))';

PHIx=KxXgrid*W*Vm*diag(Em.^(-1)); 

end

