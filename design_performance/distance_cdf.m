function [ F,x ] = distance_cdf( Xn, Xtest, k_fracfact )
% function [ F,x ] = distance_cdf( Xn, Xtest, k_fracfact )
%__________________________________________________________________________
% For Xn (d*n) an n-point design in [0,1]^d, returns the normalized
% the empirical distance cdf F(x), i.e., the cdf of n^(1/d)*d(U,Xn), where
% d(U,Xn) = min_i ||U-Xn(:,i)||, with U uniformly distributed on the set
% Xtest, a d*Q matrix, consisting for instance of scrambled Sobol' points 
% when the space is [0,1]^d:
%   pSmm = sobolset(d); pSmm = scramble(pSmm,'MatousekAffineOwen');
%   Xtest=(net(pSmm,2^19))';
% When k_fracfact=k>0, Xtest is completed by a k^d fractional factorial design 
%      (==> k should be kept small form large d, k=2 is reasonable)
% ---> one may use [ Xtest ] = candidate_set( d,q,k,center,Scramble )
% F and x are column vectors of length Q
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE.

[d,n]=size(Xn);
if isempty(Xtest)==0 && k_fracfact>0
    % complete by a k^d fractional factorial design
    k=k_fracfact;
    Complement=((fullfact(k*ones(1,d))-1)/(k-1))';
    Xtest=[Xtest Complement]; Xtest=(unique(Xtest',"rows"))';
end
Vd=1;
DD=(n*Vd)^(1/d)*pdist2(Xn',Xtest');
DXn=min(DD,[],1);
[F,x]=ecdf(DXn);

end