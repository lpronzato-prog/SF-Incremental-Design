function [ X, Gen, Dist_pattern, WLP ] = greedy_generators( d,m,G_initial,nrandom,excursion,iterations )
% function [ X, Gen, Dist_pattern, WLP ] = greedy_generators( d,m,G_initial,random,excursion,iterations )
%__________________________________________________________________________
% Greedy construction of generating equations for a 2^(d-m) fractional
% factorial design based on maximin distance
% Returns X, the list of generators in matrix form), the distance distribution 
% (pattern) and the word length pattern WLP
% G_initial (n_initial*(d-m)) contains a list of n_initial generators 
% in matrix form, 
% e.g., G_initial(1,:)=[0 0 1 2 4 5]; G_initial(2,:)=[0 2 3 4 5 6]; for d-m=6
% with all elements in {0,...,d-m} and n_initial <= m
% An example is: [X,Gen]=greedy_generators( 16,8,...
%                           [1 2 3 4 5 6 7 8
%                            0 0 0 1 2 3 4 5],0,1,10 );
%
% --> [ X, Gen, Dist_pattern, WLP ]= greedy_generators( d,m,Gen,0,0,1)
%       simply "evaluates" the design corresponding to generators in Gen
% --> [~,Gen,~,~]=greedy_generators(d,m,[],m,0,1) 
%       generates a set Gen (in matrix form) of m random generators
%
% Uses compare_patterns.m
%
% If the integer "iterations" is > 1, an algorithm of Mitchell's DETMAX type 
% [T.J. Mitchell, An algorithm for the construction of D-optimal experimental 
% designs, Technometrics, 16:203-210, 1974] is used to find good (hopefully 
% optimal) generators, with excursions going up to "excursion" beyond m
% If nrandom>0, the initialization is by a random choice of nrandom generators,
% for example: greedy_generators( 16,8,[],2,1,10 );
%-----
% Author: L. Pronzato, 2020 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE.

k=d-m;
X0=2*(fullfact(2*ones(1,k))-3/2); % full factorial 2^k design in {-1,1}^k

if excursion > 2^k-k-1-m
    disp('excursion cannot be too large')
    excursion=2^k-k-1-m
end    

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

Gen=zeros(m+excursion,k);
X_add=zeros(2^k,m+excursion);
n_initial=0;
dim=k;

nrandom=min(nrandom,m+excursion);

if nrandom>0
    iselect=randperm(n_generators,nrandom);
    G_initial=ALL_GEN(iselect,:);
end    

if isempty(G_initial)==0
    G_initial=sort(G_initial,2); % order the factors is ascending order
    [n_initial,~]=size(G_initial);
    Gen(1:n_initial,:)=G_initial;
    for j=1:n_initial
        gen_j=G_initial(j,:); gen_j=gen_j(gen_j>0);
        X_add(:,j)=prod(X0(:,gen_j),2);
    end
    dim=k+n_initial; % dim is only used to compute Hamming distance in Matlab
end  

iterations=max(iterations,1); % To run the algorithm at least once 
                              % and construct one set of m generators 
                              % in any case

for iter=1:iterations                              
    for i_add=n_initial+1:m+excursion
        dim=k+i_add;
        Best_Dist_pattern=2^k*ones(1,dim+1); % to keep trace of best distance pattern
        for j=1:n_generators
            Gen_j=ALL_GEN(j,:); % a generator
            if ismember(Gen_j,Gen,'rows')==0
                % this is really a new generator to be considered
                gen_j=Gen_j(Gen_j>0); 
                X_j=prod(X0(:,gen_j),2); % vector of coordinates in additional dimension
                Xtest=[X0 X_add(:,1:i_add-1) X_j];
                DH_j=dim*pdist2(Xtest(1,:),Xtest,'hamming');
                Dist_pattern=hist(DH_j,0:1:dim);
                ibest=compare_patterns( Dist_pattern,Best_Dist_pattern );
                if ibest==1
                    Best_Dist_pattern=Dist_pattern;
                    Gen(i_add,:)=Gen_j;
                    X_add(:,i_add)=X_j;
                end
            end
        end
    end
    % We have a "greedy" design with m+excursion generators
    if excursion>0        
        % we must return to a design with m generators
        % --> remove "excursion" generators from Gen
        % --> remove "excursion" columns to X
        Gen_removed=zeros(excursion,k); % the generators removed
        Gen_final=Gen(m+1:m+excursion,:);
        nb_removed=0;
        for i_add=m+excursion:-1:m+1
            % find worst generator
            dim=k+i_add-1;  % dimension once an additional generators will 
                            % have been removed
            Best_Dist_pattern=2^k*ones(1,dim+1); % to keep trace of best distance
            for j=1:i_add % test all generators that remain
                i_remain=setdiff(1:i_add,j);
                Xtest=[X0 X_add(:,i_remain)];
                DH_j=dim*pdist2(Xtest(1,:),Xtest,'hamming');
                Dist_pattern=hist(DH_j,0:1:dim);
                ibest=compare_patterns( Dist_pattern,Best_Dist_pattern );
                if ibest==1
                    Best_Dist_pattern=Dist_pattern;
                    j_removed=j;
                end
            end
            % keep trace of generators removed
            nb_removed=nb_removed+1;
            Gen_removed(nb_removed,:)=Gen(j_removed,:);
            % remove generator and column j_removed from Gen and X_add
            Gen(j_removed,:)=[];
            X_add(:,j_removed)=[];
        end
        n_initial=m;
        G_initial=Gen;
        if isempty(setdiff(Gen_final,Gen_removed,'rows')) && iter>= 2
            break % we remove the last generators introduced
        end 
    end
end
X=[X0 X_add(:,:)];
% THIS IS THE END...
DH=dim*pdist2(X(1,:),X,'hamming');
Dist_pattern=hist(DH,0:1:dim);
WLP=From_distancepattern_to_WLP( Dist_pattern,d,m );
end

