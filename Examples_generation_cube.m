% Examples of design generations in [0,1]^d
% The design Xi is generated with Algorithm 6.i of [Karvonen, Pronzato and Zhigljavsky, 2026]

close all
clear variables

oldpath = path; path(oldpath,'auxiliary_functions')
newpath = path; path(newpath,'design_generation')
newpath = path; path(newpath,'design_performance')

%--------------------------------------------------------------------------
d=2             % dimension
n=25            % length of the design sequence 

C=11; k=2; center=1; Xcand=candidate_set( d,C,k,center,0 ); % = candidate set, with 2^C points
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Greedy packing
[ X1,Time ] = coffee_house_adapt( Xcand, n, zeros(1,n), 1, 1);
figure(1)
plot([0 1 1 0 0],[0 0 1 1 0],'k-',X1(1,:),X1(2,:),'r.','LineWidth',2,'MarkerSize',15),  axis('square'), axis('equal'), axis off

%--------------------------------------------------------------------------
% Boundary-phobic greedy packing
    % like in [Nogalez-Gómes, Pronzato and Rendas, 2021]
    Vd=pi^(d/2)/gamma(d/2+1); 
    beta=d/(2*(n*Vd)^(-1/d))-sqrt(d);        
    Alpha=(1/beta)*ones(1,n);
    % like in [Shang and Apley 2020]
    beta=2*sqrt(2*d);
    Alpha=(1/beta)*ones(1,n);
[ X1_beta,Time ] = coffee_house_adapt( Xcand, n, Alpha, 1, 1);
figure(11)
plot([0 1 1 0 0],[0 0 1 1 0],'k-',X1_beta(1,:),X1_beta(2,:),'r.','LineWidth',2,'MarkerSize',15),  axis('square'), axis('equal'), axis off

%--------------------------------------------------------------------------
% relaxed greedy packing
    % based on [Nogalez-Gómes, Pronzato and Rendas, 2021]
    Vd=pi^(d/2)/gamma(d/2+1); Rn=(n*Vd)^(-1/d);
    Alpha=(1-2*Rn/sqrt(d))*ones(1,n);
[ X2,Time ] = coffee_house_relaxed( Xcand,n,Alpha,0 );
figure(2)
plot([0 1 1 0 0],[0 0 1 1 0],'k-',X2(1,:),X2(2,:),'r.','LineWidth',2,'MarkerSize',15),  axis('square'), axis('equal'), axis off

%--------------------------------------------------------------------------
% Greedy packing with projections
    dims=1:d; Wdim=1:d;
[ X1_PR,Time ] = coffee_house_subspaces( Xcand, n, dims, Wdim, 0, 0);
figure(111)
plot([0 1 1 0 0],[0 0 1 1 0],'k-',X1_PR(1,:),X1_PR(2,:),'r.','LineWidth',2,'MarkerSize',15),  axis('square'), axis('equal'), axis off

%--------------------------------------------------------------------------
% Greedy energy minimization
    kernel='Riesz';   
    tensorised=0;    
    corr_length=NaN; 
    power=2*d;
    tau=0;
[ X3,Time ] = kernel_herding_Energy_weighted_MMD( Xcand, n, NaN,power,kernel,tensorised,tau, [] );
figure(3)
plot([0 1 1 0 0],[0 0 1 1 0],'k-',X3(1,:),X3(2,:),'r.','LineWidth',2,'MarkerSize',15),  axis('square'), axis('equal'), axis off

%--------------------------------------------------------------------------
% Greedy MMD minimization
    kernel='matern32';
    tensorised=1;
    corr_length=n^(-1/d); 
    power=NaN;
    tau=1;
[ X4,Time ] = kernel_herding_Energy_weighted_MMD( Xcand, n, 1/corr_length,power,kernel,tensorised,tau, [] );
figure(4)
plot([0 1 1 0 0],[0 0 1 1 0],'k-',X4(1,:),X4(2,:),'r.','LineWidth',2,'MarkerSize',15),  axis('square'), axis('equal'), axis off

%--------------------------------------------------------------------------
% Greedy minimization of the double relaxation of covering
    kernel='riesz';
    tensorised=0;    
    corr_length=NaN; 
    power=2*d;
% simple greedy
    %[ X5,Time ] = greedy_relaxed_CRKmu( Xcand, n, 1/corr_length,power,kernel,tensorised,Xcand );
% lazy-greedy version
    [ X5,Time ] = lazy_greedy_relaxedCRKmu( Xcand, n, 1/corr_length,power,kernel,tensorised,Xcand );
figure(5)
plot([0 1 1 0 0],[0 0 1 1 0],'k-',X5(1,:),X5(2,:),'r.','LineWidth',2,'MarkerSize',15),  axis('square'), axis('equal'), axis off
    
%--------------------------------------------------------------------------
% Greedy minimization of the Ls-mean quantization error   
    s=2*d;
% simple greedy
    %[ X6,Time ] = greedy_quantization( Xcand, Xcand, n, s);    
% lazy-greedy version
    [ X6,Time ] = lazy_greedy_quantization( Xcand, Xcand, n, s);
figure(6)
plot([0 1 1 0 0],[0 0 1 1 0],'k-',X6(1,:),X6(2,:),'r.','LineWidth',2,'MarkerSize',15),  axis('square'), axis('equal'), axis off

%--------------------------------------------------------------------------
% Maximum Entropy Sampling        
    kernel='matern32';  
    tensorised=1;    
    corr_length=n^(-1/d); 
    power=NaN;
[ X7,Time ] = greedy_MES( Xcand, n, 1/corr_length,power,kernel,tensorised);    
figure(7)
plot([0 1 1 0 0],[0 0 1 1 0],'k-',X7(1,:),X7(2,:),'r.','LineWidth',2,'MarkerSize',15),  axis('square'), axis('equal'), axis off

%--------------------------------------------------------------------------
% Greedy minimization of the IMSPE
    kernel='matern32';   
    tensorised=1;    
    corr_length=n^(-1/d); 
    power=NaN;
[ X8,Time ] = greedy_IMSPE( Xcand, n, 1/corr_length,power,kernel,tensorised,[]);
figure(8)
plot([0 1 1 0 0],[0 0 1 1 0],'k-',X8(1,:),X8(2,:),'r.','LineWidth',2,'MarkerSize',15),  axis('square'), axis('equal'), axis off

%--------------------------------------------------------------------------
% Greedy minimization of the IMSPE with spectral decomposition
    kernel='matern32';   
    tensorised=1;    
    corr_length=n^(-1/d); 
    power=NaN;
    mm=3*n; n1=100; m1=100; w1=ones(1,n1)/n1;
[ X9,Time ] = greedy_IMSPE_KarhunenL( Xcand, n, mm,m1,n1,w1,1/corr_length,power,kernel);
figure(9)
plot([0 1 1 0 0],[0 0 1 1 0],'k-',X9(1,:),X9(2,:),'r.','LineWidth',2,'MarkerSize',15),  axis('square'), axis('equal'), axis off

%--------------------------------------------------------------------------
% Greedy minimization of the IMSPE based on a Bayesian Linear Model;
    kernel='matern32';   
    tensorised=1;    
    corr_length=n^(-1/d); 
    power=NaN;
    mm=3*n; n1=100; m1=100; w1=ones(1,n1)/n1;
[ X10,Time ] = greedy_IMSPE_BLM( Xcand, n, mm,m1,n1,w1,1/corr_length,power,kernel);
figure(10)
plot([0 1 1 0 0],[0 0 1 1 0],'k-',X10(1,:),X10(2,:),'r.','LineWidth',2,'MarkerSize',15),  axis('square'), axis('equal'), axis off
