function [ alphaCR ] = Covering_quantile( Xn, sequence, Xtest, k_fracfact, alpha )
% function [ alphaCR ] = Covering_quantile( Xn, sequence, Xtest, k_fracfact, alpha )
%__________________________________________________________________________
% For Xn (d*n) an n-point design, returns 
% alphaCR of size m = length(sequence), a row vector containing
% the alpha-covering quantiles of the distance of random point in in [0,1]^d
% to the sequence of designs Xn(:,(1:i1)),...,Xn(:,(1:im)), with 
% [i_1 ... i_m] = sequence = a vector of m consecutive integers in {1,...,n}
% If sequence = k, alphaCR is just the alpha-covering quantile of Xn(:,(1:k))
% The alpha-covering quantiles are estimated by evaluation on the Q points in Xtest, 
% a d*Q matrix, for instance of scrambled Sobol' points:
%   pSmm = sobolset(d); pSmm = scramble(pSmm,'MatousekAffineOwen');
%   Xtest=(net(pS,2^19))';
% When k_fracfact=k>0, Xtest is completed by a k^d fractional factorial design 
%      (==> k should be kept small form large d, k=2 is reasonable)
% ---> one may use [ Xtest ] = candidate_set( d,q,k,center,Scramble )
%-----
% Author: L. Pronzato, 2021 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

nmax=sequence(end);
[d,n]=size(Xn);
if nmax>n
    display('sequence is too long for the nb. of design points Xn') 
    alphaCR=NaN(size(sequence));
    return
end

if isempty(Xtest)==0 && k_fracfact>0
    % complete by a k^d fractional factorial design
    k=k_fracfact;
    Complement=((fullfact(k*ones(1,d))-1)/(k-1))';
    Xtest=[Xtest Complement];
end
[~,N]=size(Xtest);
index_left_alpha_quantile=max(ceil(alpha*N),1);
LalphaCR=length(sequence);
alphaCR=NaN(1,LalphaCR);
timebar= waitbar(0,'covering quantile...'); 
Dist_all=pdist2(Xn(:,1:nmax)',Xtest');
for i=1:LalphaCR
    waitbar(i/LalphaCR,timebar);
    ialphaCR=sequence(i);
    if i==1
        D2Xn=min(Dist_all(1:ialphaCR,:),[],1);
    else
        D2Xn=min(D2Xn,Dist_all(ialphaCR,:));
    end    
    D2Xnsort=sort(D2Xn); 
    alphaCR(i)=D2Xnsort(index_left_alpha_quantile);
end
close(timebar)
end
