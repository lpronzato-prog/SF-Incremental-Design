function [ discr ] = discrepancy_wa( X )
% function [ discr ] = discrepancy_wa( X )
%__________________________________________________________________________
% For X (d*n) an n-point design in [0,1]^d, returns 
% the wrap-around discrepancy [Fang, Li, Sudjianto, 2006, page 73]
% = MMD for the tensor-product kernel K(x,x')=3/2-|x-x'|(1-|x-x'|)
% for X a d*n design matrix made of n points in [0,1]^d
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE.

[d,n]=size(X);
D=zeros(n,n);
for i=1:n
    for j=i:n
        Dij=abs(X(:,i)-X(:,j));
        D(i,j)=prod(1.5*ones(d,1)-Dij.*(ones(d,1)-Dij));
	end
end
D=D+triu(D,1)';
T2=sum(sum(D))/n^2;
discr=sqrt(-(4/3)^d+T2);
end

