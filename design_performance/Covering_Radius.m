function [ CR ] = Covering_Radius( Xn, sequence, d_threshold, Xtest, k_fracfact )
% function [ CR ] = Covering_Radius( Xn, sequence, d_threshold, Xtest, k_fracfact )
%__________________________________________________________________________
% For Xn (d*n) an n-point design, returns 
% CR of size m = length(sequence), a row vector containing the covering 
% radii of the sequence of designs Xn(:,(1:i1)),...,Xn(:,(1:im)), with 
% [i_1 ... i_m] = sequence = a vector of m consecutive integers in {1,...,n}
% If sequence = k, CR is just CR(Xn(:,(1:k)))
%
% - for d<=d_threshold, the exact value is calculated (via Voronoi diagram)
%       for the design space XX=[0,1]^d -- uses voronoi_truncated01.m
% - for d>d_threshold, CR is (under-)estimated by evaluation on the Q points 
%       in Xtest --- which can thus be used for any design space
% Xtest is a d*Q matrix, for instance of scrambled Sobol' points when XX=[0,1]^d:
%   pSmm = sobolset(d); pSmm = scramble(pSmm,'MatousekAffineOwen');
%   Xtest=(net(pS,2^19))';
% To keep the computational cost reasonable, d_threshold should be <= 5 (this is enforced)
% When k_fracfact=k>0, Xtest is completed by a k^d fractional factorial design 
%      (==> k should be kept small form large d, k=2 is reasonable)
%-----
% Author: L. Pronzato, 2021 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

[d,~]=size(Xn);
d_threshold=min(d_threshold,5);
%
nmax=sequence(end);
[~,n]=size(Xn);
if nmax>n
    display('sequence is too long for the nb. of design points Xn') 
    CR=NaN(size(sequence));
    return
end

if isempty(Xtest)==0 && k_fracfact>0
    % complete by a k^d fractional factorial design
    k=k_fracfact;
    Complement=((fullfact(k*ones(1,d))-1)/(k-1))';
    Xtest=[Xtest Complement];
end
LCR=length(sequence);
CR=NaN(1,LCR);
timebar= waitbar(0,'covering radius...'); 
if d>d_threshold
    % Compute all what we need once for all
    Dist_all=pdist2(Xn(:,1:nmax)',Xtest');
end    
for i=1:LCR
    waitbar(i/LCR,timebar);
    iCR=sequence(i);
    if d==1
        CR(i) = Minimax_d1( Xn(:,1:iCR), 0,1 );
    elseif d<=d_threshold
        [ CR(i), ~, ~, ~ ] = voronoi_truncated01( Xn(:,1:iCR), [], [], 0 );
    else
        if i==1
            D2Xn=min(Dist_all(1:iCR,:),[],1);
        else
            D2Xn=min(D2Xn,Dist_all(iCR,:));
        end    
        CR(i)=max(D2Xn);
    end
end
close(timebar)
end
