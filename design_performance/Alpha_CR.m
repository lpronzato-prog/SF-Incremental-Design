function [ Alpha,CR ] = Alpha_CR( Xn, d_threshold, Xtest, k_fracfact )
% function [ Alpha,CR ] = Alpha_CR( Xn, d_threshold, Xtest, k_fracfact )
%__________________________________________________________________________
% For Xn (d*n) an n-point design, returns 
% a (1*n) row vector CR containing the covering radii of the designs Xn(:,(1:i)), i=1,...,n, 
% computed by CR = Covering_Radius( Xn, 1:n, d_threshold, Xtest, k_fracfact )
% and PR and CR efficiency bounds: the PR and CR efficiencies (evaluated on Xtest) 
% of all Xn(:,(1:k)) are at least ( min_{i=1,...,k} Alpha(i) )/2
% where Alpha(i)=min_{x in Xn(:,(1:i))} = ||Xn(:,(1:i+1)-x||/CR[Xn(:,(1:i))]
% is the ith component of the 1*(n-1) vector Alpha
%
% for d<=d_threshold, the exact CR value is calculated (via Voronoi diagram)
% for d>d_threshold, CR is underestimated by evaluation on the Q points in Xtest, 
% a d*Q matrix, for instance of scrambled Sobol' points when the space is [0,1]^d:
%   pSmm = sobolset(d); pSmm = scramble(pSmm,'MatousekAffineOwen');
%   Xtest=(net(pS,2^19))';
% To keep the computational cost reasonable, d_threshold should be <= 5 (this is enforced)
% When k_fracfact=k>0, Xtest is completed by a k^d fractional factorial design 
%      (==> k should be kept small form large d, k=2 is reasonable)
% ---> one may use [ Xtest ] = candidate_set( d,q,k,center,Scramble )
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

[d,n]=size(Xn);
CR = Covering_Radius( Xn, 1:n, d_threshold, Xtest, k_fracfact );
Dnp1=min(squareform(pdist(Xn'))+sqrt(d)*triu(ones(n,n)),[],2); 
Alpha=Dnp1(2:n)'./CR(1:n-1);
end
