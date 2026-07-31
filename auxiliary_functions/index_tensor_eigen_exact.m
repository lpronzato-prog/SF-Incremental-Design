function [ MM,T, Lambda ] = index_tensor_eigen_exact( lambda,d, nTmax, degree_max )
% function [ MM,T, Lambda ] = index_tensor_eigen_exact( lambda,d, nTmax, degree_max )
%__________________________________________________________________________  
% lambda a vector m*1 sorted by decreasing values, 
% requires : d>=2, m^d>=nTmax
% returns a matrix MM (L*m) of integers with row i of the form
% MM(i,:)=(k_1, ..., k_m) and such that 
% Lambda(i)=lambda.^(k_1, ..., k_m)=prod_j=1^m lambda(j)^k_j<=Lambda(i+1)
% and sum_j=1^m k_j=d
% The corresponding component T(i) of the vector T gives the number of
% combinations that will have eigenvalue Lambda(i)
% 
% 1/ There are T(i)=d!/(k_1! k_2! ... k_m!) possible tensorisations
% corresponding to MM(i,:) (they are obtained by tensorised_row.m)
% 2/ L is the smallest integer such that SumTi=sum_{i=1}^L T_i >= nTmax. 
% 3/ The Lambda(i) are decreasing
%
% Set degree_max = Inf to base the construction on the selection of nTmax
% largest eigenvalues only
%-----
% L. Pronzato, 2025
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
%-----

lambda=lambda(:); 
m=length(lambda);
factd=factorial(d);
nTmax=min(nTmax,m^d);

% Initialisation  
MM=[d zeros(1,m-1)]; 
T=1; % T(i)=T_i for row i of MM
Lambda=lambda(1)^d;
SumTi=1; 

MM_descendants=[d-1 1 zeros(1,m-2)]; % The initial MM has one descendant only 
Lambda_descendants=lambda(1)^(d-1)*lambda(2);
T_descendants=d; %( = d!/(d-1)!)

while SumTi<nTmax
    % add descendant with largest Lambda to MM
    [~,imax]=max(Lambda_descendants);
    MM_add=MM_descendants(imax,:);
    if MM_add*(0:1:m-1)'<=degree_max % Used for polynomial chaos models
        MM=[MM;MM_add];
        % update also T and Lambda
        T=[T; T_descendants(imax)];
        Lambda=[Lambda;Lambda_descendants(imax)];
    end
    SumTi=SumTi+T_descendants(imax);
    
    % Remove the corresponding rows from all _descendants
    MM_descendants(imax,:)=[];
    T_descendants(imax)=[];
    Lambda_descendants(imax)=[];
    
    if SumTi<nTmax % useless to do the job otherwise
        % Compute descendants of MM_add
        iopen=find(MM_add~=0); % indices j such that k_j > 0
        for j=1:length(iopen)
            jj=iopen(j); % MM_add(k,jj)>0
            if jj<m 
                % try to add a row to MM_add by shifting one eigenvalue to the
                % right
                MMj=MM_add; MMj(jj)=MMj(jj)-1; MMj(jj+1)=MMj(jj+1)+1;
                if ismember(MMj,MM_descendants,'rows')==0 % The row may already be present
                    MM_descendants=[MM_descendants
                                    MMj];
                    % Also add a row to T_descendants and Lambda_descendants
                    LTj=log(factd)-sum(log(factorial(MMj(1:m)))); Tj=round(exp(LTj));
                    T_descendants=[T_descendants;
                                        Tj];
                    Lambdaj=exp(MMj*log(lambda));
                    Lambda_descendants=[Lambda_descendants;
                                        Lambdaj];
                end
            end
        end
    end
end

end

