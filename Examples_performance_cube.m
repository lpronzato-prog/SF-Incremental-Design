% Examples of design performance evaluation in [0,1]^d
% for some of the criteria considered in [Karvonen, Pronzato and Zhigljavsky, 2026]

close all
clear variables

oldpath = path; path(oldpath,'auxiliary_functions')
newpath = path; path(newpath,'design_generation')
newpath = path; path(newpath,'design_performance')

%--------------------------------------------------------------------------
d=2             % dimension
n=25            % design size

Xn=rand(d,n);    % a random design

% variables used for performance evaluation
    C=15; k=2; center=1; Xtest=candidate_set( d,C,k,center,1 ); % = candidate set, with 2^C points
    kernel='matern32';
    corr_length=n^(-1/d); 
    power=NaN;
    tensorised=1;
%--------------------------------------------------------------------------

% Covering radius
%   [ CR ] = Covering_Radius( Xn, 1:n, 3, Xtest, 0 ) 
% Packing radius
%   [ PR ] = Packing_Radius( Xn, 1:n )
% PR and CR efficiency bounds
%   [ Alpha,CR ] = Alpha_CR( Xn, 3, Xtest, 0 ) 
% Ls-mean quantization error
%   [ LsME ] = Ls_Mean_Error( Xn, 1:n, Xtest, 0, 2 ) 
% alpha-covering quantile
%   [ alphaCR ] = Covering_quantile( Xn, 1:n, Xtest, 0, 0.999 ) 
% L2 symmetric discrepancy
%   [ DL2sym ] = DL2sym_1n( Xn )
% Energy and physical energy
%   [ Energy, Physical_energy ] = Energy_Xn( Xn, 1/corr_length,power,kernel,tensorised)
% Squared MMD
% [ mmd2 ] = MMD2( Xn, 1:n, 1/corr_length,power,kernel) % squared MMD
% IMSPE (exact)
    [ IMSPE ] = IntegratedMSPE( Xn, 1:n, 1/corr_length,power,kernel) 
% IMSPE (empirical)    
    [ IMSPE_2 ] = IntegratedMSPE_discrete_mu( Xn, 1:n, 1/corr_length,power,kernel,tensorised, Xtest, 0 )
% IMSPE (based on a spectral decomposition)    
    mm=3*n; n1=100; m1=100; w1=ones(1,n1)/n1; 
    [ IMSPE_3 ] = IntegratedMSPE_KarhunenL( Xn, 1:n, mm,m1,n1,w1,1/corr_length,power,kernel)
% IMSPE (based on a Bayesian Linear Model)  
    [ IMSPE_4 ] = IntegratedMSPE_BLM( Xn, 1:n, mm,m1,n1,w1,1/corr_length,power,kernel)
% Maximum of squared prediction error
%   MMSPE = MaxMSPE( Xn, 1:n, Xtest, 0, 1/corr_length,power,kernel,tensorised)
% Determinant of kernel matrix
%   detK = DetKn( Xn, 1:n, 1/corr_length,power,kernel,tensorised)
% Lebesgue constant
%   Lct = Lebesgue_ct( Xn, 1:n, Xtest, 0, 1/corr_length,power,kernel,tensorised)
% (normalized-) distance cdf
%   [ F,x ] = distance_cdf( Xn, Xtest, 0 );
% Star and extreme discrepancies of the max. distance to the center of the cube
%   [ discr_star_n, discr_n ] = star_extreme_radial_discr_sequence( Xn )

