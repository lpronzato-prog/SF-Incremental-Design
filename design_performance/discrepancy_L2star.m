function [ discr ] = discrepancy_L2star( X )
% function [ discr ] = discrepancy_L2star( X )
%__________________________________________________________________________
% For X (d*n) an n-point design in [0,1]^d, returns 
% the L2-star discrepancy (anchored at the origin),
% = MMD for the tensor-product kernel K(x,x')=1-max(x,x')
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
        D(i,j)=prod( ones(d,1)-max(X(:,i),X(:,j)) );
	end
end
D=D+triu(D,1)';
T3=sum(sum(D))/n^2;
T2=-(2/n)*sum( prod((ones(d,n)-X.^2)/2) );
discr=sqrt((1/3)^d+T2+T3);
end

