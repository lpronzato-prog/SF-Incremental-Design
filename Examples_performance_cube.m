% Examples of performance evaluation of a given design in [0,1]^d
% for some of the criteria considered in the book [Karvonen, Pronzato and Zhigljavsky, Springer, 2026]
% (the Section where the performance criterion is defined is indicated)

close all
clear variables

oldpath = path; path(oldpath,'auxiliary_functions')
newpath = path; path(newpath,'design_generation')
newpath = path; path(newpath,'design_performance')

%--------------------------------------------------------------------------
d=2             % dimension
n=25            % design size
Xn=rand(d,n);    % a random design
%--------------------------------------------------------------------------
% variables used for some of the performance criteria below
    C=15; k=2; center=1; Xtest=candidate_set( d,C,k,center,1 ); % = a test set, with 2^C points
    kernel='matern32';
    corr_length=n^(-1/d); 
    power=NaN;
    tensorised=1;
%--------------------------------------------------------------------------
% Covering radius: --> Section 1.1
   [ CR ] = Covering_Radius( Xn, 1:n, 3, Xtest, 0 ) 
% Packing radius: --> Section 1.1
   [ PR ] = Packing_Radius( Xn, 1:n )
% PR and CR efficiency bounds (given by min(Alpha)/2, --> Corollary 6.9): 
   [ Alpha,CR ] = Alpha_CR( Xn, 3, Xtest, 0 ) 
% (normalized-) distance cdf: --> Section 1.2
   [ F,x ] = distance_cdf( Xn, Xtest, 0 );
% gamma-covering quantile: --> Section 1.2.2
   [ gammaCR ] = Covering_quantile( Xn, 1:n, Xtest, 0, 0.999 ) 
% Ls-mean quantization error: --> Section 1.3
   [ LsME ] = Ls_Mean_Error( Xn, 1:n, Xtest, 0, 2 ) 
%---------------------------------   
% Energy and physical energy: Section --> 3.1.3
   [ Energy, Physical_energy ] = Energy_Xn( Xn, 1/corr_length,power,kernel,tensorised)  
% Squared MMD: Section --> 3.4.3 
 [ mmd2 ] = MMD2( Xn, 1:n, 1/corr_length,power,kernel)   
% L2 symmetric discrepancy (for all X1,...,Xn): --> Section 3.4.7
   [ DL2sym ] = DL2sym_1n( Xn )
% L2 centered discrepancy: --> Section 3.4.7   
   [ discr_c ] = discrepancy_cent( Xn )
% L2 extreme discrepancy: --> Section 3.4.7   
   [ discr_e ] = discrepancy_L2extreme( Xn )
% L2 star discrepancy (anchored at the origin): --> Section 3.4.7   
   [ discr_s ] = discrepancy_L2star( Xn )
% L2 star discrepancy (modified): --> Section 3.4.7   
   [ discr_sm ] = discrepancy_L2starM( Xn )
% L2 wrap-around discrepancy: --> Section 3.4.7      
   [ discr_wa ] = discrepancy_wa( Xn )
%---------------------------------
% Lebesgue constant: --> Section 5.2.4
   Lct = Lebesgue_ct( Xn, 1:n, Xtest, 0, 1/corr_length,power,kernel,tensorised)  
% Maximum of squared prediction error: --> Section 5.6.2 
   MMSPE = MaxMSPE( Xn, 1:n, Xtest, 0, 1/corr_length,power,kernel,tensorised)
% Determinant of kernel matrix: --> Section 5.6.2 
   detK = DetKn( Xn, 1:n, 1/corr_length,power,kernel,tensorised)  
% IMSPE (exact): --> Section 5.7.2 
    [ IMSPE ] = IntegratedMSPE( Xn, 1:n, 1/corr_length,power,kernel) 
% IMSPE (empirical): --> Section 5.7.2     
    [ IMSPE_2 ] = IntegratedMSPE_discrete_mu( Xn, 1:n, 1/corr_length,power,kernel,tensorised, Xtest, 0 )
% IMSPE (based on a spectral decomposition): --> Section 5.7.2     
    mm=3*n; n1=100; m1=100; w1=ones(1,n1)/n1; 
    [ IMSPE_3 ] = IntegratedMSPE_KarhunenL( Xn, 1:n, mm,m1,n1,w1,1/corr_length,power,kernel)
% IMSPE (based on a Bayesian Linear Model): --> Section 5.7.2   
    [ IMSPE_4 ] = IntegratedMSPE_BLM( Xn, 1:n, mm,m1,n1,w1,1/corr_length,power,kernel)
%---------------------------------    
% Star and extreme discrepancies of the max. distance to the center of the cube: --> Section 6.2.6
   [ discr_star_n, discr_n ] = star_extreme_radial_discr_sequence( Xn )

