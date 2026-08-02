function [ X, Gen, Dist_pattern, WLP ] = SA_exchange_generators( d,m,G_initial,iterations,TYPE,a, stop )
% function [ X, Gen, Dist_pattern, WLP ] = SA_exchange_generators( d,m,G_initial,iterations,TYPE,a, stop )
%__________________________________________________________________________
% A simulated annealing algorithm of the exchange type for the
% construction of generating equations for a 2^(d-m) fractional
% factorial design based on maximin distance
% = Algorithm 7.2 of [Karvonen, Pronzato and Zhigljavsky, 2026]
% Returns X, the list of generators in matrix form), the distance distribution 
% (pattern) and the word length pattern WLP
% G_initial (m*(d-m)) contains a list of m generators in matrix form, 
% --> use [~,Gen,~,~]=greedy_generators(d,m,[],m,0,1) 
%       to generate a set Gen (in matrix form) of m random generators
%
% Uses From_distancepattern_to_WLP.m, compare_patterns.m
%
% iterations = nb. of SA iterations
% TYPE=1 : -> optimizes the distance pattern
% TYPE=2 : -> optimizes the word length pattern WLP
% a --> type of temperature decrease in SA: 
%               T=T0/iter^a if a>0,
%               T=T0/log(iter+1) if a=0
% if stop=1, the algorithm stops as soon as a design with minimum hamming
%           distance >= d/4 is found (see [Cabral-Farias, Pronzato,
%           Rendas,2023] for a justification)
%-----
% Author: L. Pronzato, 2020 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE.

m_max=floor(d-log2(d+1));
if m>m_max
    'm is too large'
    X=[]; Gen=[]; Dist_pattern=[]; WLP=[];
    return
end    

k=d-m;
X0=2*(fullfact(2*ones(1,k))-3/2); % full factorial 2^k design in {-1,1}^k

% Create a big (2^k-k-1)*k matrix ALL_GEN containing all possible generators
n_generators=2^k-k-1;
ALL_GEN=zeros(n_generators,k);
istart=0; 
for n_factors=2:k
    n_i=nchoosek(k,n_factors);
    Li=[zeros(n_i,k-n_factors) nchoosek(1:k,n_factors)];
    ALL_GEN(istart+1:istart+n_i,:)=Li;
    istart=istart+n_i;
end 

% Evaluate the initial set of m generators
X_add=zeros(2^k,m);
G_initial=sort(G_initial,2); % order the factors is ascending order
Gen=G_initial;
for j=1:m
    gen_j=G_initial(j,:); gen_j=gen_j(gen_j>0);
    X_add(:,j)=prod(X0(:,gen_j),2);
end
dim=d; % dim is only used to compute Hamming distance in Matlab
Xtest=[X0 X_add];
DH_j=dim*pdist2(Xtest(1,:),Xtest,'hamming');
Dist_pattern=hist(DH_j,0:1:dim);
if TYPE==2
    % consider WLP instead of distance pattern
    Dist_pattern=From_distancepattern_to_WLP( Dist_pattern,d,m );
end
% keep trace of very best design and generators
very_best_Dist_pattern=Dist_pattern;
VB=very_best_Dist_pattern(2:d+1); jmin=find(VB>0,1); %disp(jmin)
very_best_generators=G_initial;
very_best_X_add=X_add;

if stop==1 && jmin>=d/4
    % We stop here
    iterations=0;
end    

T0=2^k;
Decrease=zeros(iterations,1);
Proba=zeros(iterations,1);
Temp=zeros(iterations,1);
timebar= waitbar(0,'SA Computing...'); 
for iter=1:iterations
    if rem(iter,1000)==0
        waitbar(iter/iterations,timebar)
    end
    if a==0
        T=T0/log(iter+1);
    else
        T=T0/iter^a;
    end
    Temp(iter)=T;
    % choose a generator at random among current ones
    i_exchange=randperm(m,1);
    % choose another one at random
    repeat=1;
    while repeat==1
        j_exchange=randperm(n_generators,1);
        % make sure j_exchange does not correspond to a generator in Gen
        Gen_j_exchange=ALL_GEN(j_exchange,:);
        if ismember(Gen_j_exchange,Gen,'rows')==0
            repeat=0;
        end
    end
    % try the exchange of i_exchange by j_exchange
    gen_j_exchange=Gen_j_exchange(Gen_j_exchange>0); 
    X_j_exchange=prod(X0(:,gen_j_exchange),2); % vector of coordinates for j_exchange
    X_add_exchange=X_add;
    X_add_exchange(:,i_exchange)=X_j_exchange;
    Xtest=[X0 X_add_exchange];
    DH_exchange=dim*pdist2(Xtest(1,:),Xtest,'hamming');
    Dist_pattern_exchange=hist(DH_exchange,0:1:dim);
    if TYPE==2
        % consider WLP instead of distance pattern
        Dist_pattern_exchange=From_distancepattern_to_WLP( Dist_pattern_exchange,d,m );
    end

    % keep trace of very best design and generators
    ibest=compare_patterns( Dist_pattern_exchange,very_best_Dist_pattern );
    Gen_exchange=Gen; Gen_exchange(i_exchange,:)=ALL_GEN(j_exchange,:);
    if ibest==1
        very_best_Dist_pattern=Dist_pattern_exchange;
        VB=very_best_Dist_pattern(2:d+1); jmin=find(VB>0,1); 
        very_best_generators=Gen_exchange;
        very_best_X_add=X_add_exchange;
        
        if stop==1 && jmin>=d/4
            % We stop here
            break
        end
        
    end
    % compute the probability of acceptance
    Diff=Dist_pattern_exchange-Dist_pattern;
    kstar=find((Diff~=0),1);
    Delta=0;
    if isempty(kstar)==0
        Delta=(Dist_pattern_exchange(kstar)-Dist_pattern(kstar));
    end
    Prob=min(exp(-Delta/T),1);
    % accept the move or not...
    z=rand(1,1);
    Decrease(iter)=0;
    Proba(iter)=Prob;
    if z<=Prob
        % accept exchange
        X_add=X_add_exchange;
        Dist_pattern=Dist_pattern_exchange;
        Gen=Gen_exchange;
        Decrease(iter)=Delta;
    end
end
close(timebar)

% THIS IS THE END...

    figure(1)
    subplot(311), plot(1:iterations,cumsum(Decrease),'r-','LineWidth',2)
        title('Cumulative improvement','FontSize',20), axis tight, set(gca,'FontSize',15), grid on
    subplot(312), plot(1:iterations,Proba,'r-','LineWidth',2)
        title('Probability of acceptance','FontSize',20), axis tight, set(gca,'FontSize',15), grid on
    subplot(313), plot(1:iterations,Temp,'r-','LineWidth',2)
        title('Temperature','FontSize',20), axis tight, set(gca,'FontSize',15), grid on

Gen=very_best_generators;
X=[X0 very_best_X_add(:,:)];
Dist_pattern=very_best_Dist_pattern;

if TYPE==1
    % Compute WLP from distance pattern
    WLP=From_distancepattern_to_WLP( Dist_pattern,d,m );
else
    WLP=Dist_pattern;
    % Compute distance pattern
    DH=dim*pdist2(X(1,:),X,'hamming');
    Dist_pattern=hist(DH,0:1:dim);
end
end
