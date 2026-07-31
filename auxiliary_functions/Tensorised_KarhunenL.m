function [ LAMBDA,PHI ] = Tensorised_KarhunenL( Xcand,m,m1,n1,w1,theta,p,kernel )
% function [ LAMBDA,PHI ] = Tensorised_KarhunenL( Xcand,m,m1,n1,w1,theta,p,kernel )
%__________________________________________________________________________  
% --> uses KarhunenL_dim1.m, KarhunenL_dim1_interpolate.m,
%       index_tensor_eigen_exact.m and tensorised_row.m
%__________________________________________________________________________
% Xcand a d*n design matrix, of n design points in [0,1]^d,
% computes the first m<=m1^d eigenpairs Lambda_i and Phi_i(X)
% for a Karhunen-Loève decomposition of a tensor-product kernel and product measure.
% The measure and one-dimensional kernel (for instance kernel='Matern32') is the same 
% in all coordinates. 
% -- First, m1<=n1 eigenpairs lambda_i, phi_i are computed 
%   for this one-dimensional kernel, i=1,...,m1, using the n1 quadrature points 
%   (0:1:n1-1)/(n1-1) in [0,1] weighted by w1 (1*n1, with sum_i w1_i=1). 
% -- Canonical extensions of the phi_i are then computed for the 
%   points Xcand(k,:), for each coordinate k=1,...,d. The lambda_i and phi_i
%   are aggregated to select the first m dominant eigenpairs.
%
% In fact, we compute mp>m eigenpairs, because some eigenvalues are
%       equal: mp is such that the next one is strictly smaller 
% PHI (n*mp) contains the values of the mp eigenfunctions at the n design points
% LAMBDA (1*mp) contains the mp eigenvalues
%-----
% L. Pronzato, 2025
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.


[d,n]=size(Xcand);
if m1>n1
    display('m1 cannot be larger than n1') 
    LAMBDA=NaN; PHI=NaN;
    return
end
if m>m1^d
    display('m cannot be larger than m1^d') 
    LAMBDA=NaN; PHI=NaN;
    return
end
% Construct eigenfunctions on [0,1]
[ Vm, Em, ~ ] = KarhunenL_dim1( -1,m1,n1,w1,theta,p,kernel );
% Compute canonical extensions 
BigPhix=zeros(n,m1,d); % Contains canonical extensions computed at coordinates of Xcand
for i=1:d
    Phix_i=KarhunenL_dim1_interpolate( Vm,Em,w1,theta,p,kernel,Xcand(i,:));
    BigPhix(:,:,i)=Phix_i;
end
% Construct the tensorised form
[ MM,T,Lambda ] = index_tensor_eigen_exact( Em,d, m, Inf ); 
[LMM,~]=size(MM);
mm=sum(T);
PHI=zeros(n,mm);
LAMBDA=zeros(1,mm);
LAMBDAm1=zeros(1,mm);
nb_eigen=0;
INDICES=zeros(mm,d);    % each row (of length d) corresponds to one term of the tensorised model, 
                        % it indicates which one-dimensional (eigen-) functions are used for each variable
for j=1:LMM
    [ MTj ] = tensorised_row( MM(j,:) ); % indicates all possible combinations that can be constructed 
                                         % from a_1,...a_m, when a_i can be chosen k_i times and
                                         % MM(j,:) = [k_1,k_2,...,k_m]
    [LMTj,~]=size(MTj);
    for k=1:LMTj
        nb_eigen=nb_eigen+1; % nb of eigenpairs considered
        PHI(:,nb_eigen)=BigPhix(:,MTj(k,1),1);
        for dim=2:d    
            PHI(:,nb_eigen)=PHI(:,nb_eigen).*BigPhix(:,MTj(k,dim),dim);
        end
        INDICES(nb_eigen,:)=MTj(k,:);
        LAMBDA(nb_eigen)=Lambda(j);
        LAMBDAm1(nb_eigen)=1/Lambda(j);
    end            
end
end

