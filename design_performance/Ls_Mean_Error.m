function [ LsME ] = Ls_Mean_Error( Xn, sequence, Xtest, k_fracfact, s )
% function [ LsME ] = Ls_Mean_Error( Xn, sequence, Xtest, k_fracfact, s  )
%__________________________________________________________________________
% For Xn (d*n) an n-point design, returns 
% LsME of size m = length(sequence), a row vector containing
% the Ls-mean quantization errors,
% i.e., the Ls-norms of the distance function  
% (int_X d^s(x,X_i) dmu(x))^(1/s) with mu uniform on the design space XX,
% for the sequence of designs Xn(:,(1:i1)),...,Xn(:,(1:im)), with 
% [i_1 ... i_m] = sequence = a vector of m consecutive integers in {1,...,n}
% If sequence = k, LsME is just the Ls-mean quantization error of Xn(:,(1:k))
% The Ls-errors are estimated by evaluation on the Q grid points in Xtest,  
%   a d*Q matrix, consisting for instance when XX=[0,1]^d of scrambled Sobol' points:
%   pSmm = sobolset(d); pSmm = scramble(pSmm,'MatousekAffineOwen');
%   Xtest=(net(pSmm,2^19))';
% When k_fracfact=k>0, Xtest is completed by a k^d fractional factorial design 
%      (==> k should be kept small form large d, k=2 is reasonable)
% ---> one may use [ Xtest ] = candidate_set( d,q,k,center,Scramble )
%-----
% Author: L. Pronzato, 2022 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE.

if s<=0
    display('s must be >0')
    LsME=NaN(size(sequence));
    return
end
nmax=sequence(end);
[d,n]=size(Xn);
if nmax>n
    display('sequence is too long for the nb. of design points Xn') 
    LsME=NaN(size(sequence));
    return
end

if isempty(Xtest)==0 && k_fracfact>0
    % complete by a k^d fractional factorial design
    k=k_fracfact;
    Complement=((fullfact(k*ones(1,d))-1)/(k-1))';
    Xtest=[Xtest Complement];
end
LLsME=length(sequence);
LsME=NaN(1,LLsME);
timebar= waitbar(0,'Ls-mean quantization error...');
% Compute all what we need once for all
Dist_all_s=pdist2(Xn(:,1:nmax)',Xtest'); 
if s~=1
    Dist_all_s=Dist_all_s.^s;
end    
for i=1:LLsME
    waitbar(i/LLsME,timebar);
    i_sequence=sequence(i);
    % recursive updating of (distance to the design)^s
    if i==1
        D2Xn_s=min(Dist_all_s(1:i_sequence,:),[],1);
    else 
        D2Xn_s=min(D2Xn_s,Dist_all_s(i_sequence,:));
    end    
    LsME(i)=mean(D2Xn_s);
end
LsME=LsME.^(1/s);
close(timebar)
end
