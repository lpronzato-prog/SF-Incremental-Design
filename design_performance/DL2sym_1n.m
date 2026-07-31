function [ DL2sym ] = DL2sym_1n( X )
% function [ DL2sym ] = DL2sym_1n( X )
%__________________________________________________________________________
% For X (d*n) an n-point design in [0,1]^d, returns DL2sym, a 1*n vector of
% L2 symmetric discrepancies 
% (= MMD for the tensor-product kernel K(x,x')=2(1-|x-x'|))
% for all X(:,1:k), k=1,...,n 
% Coincides with sqrt(MMD2( Xn, sequence, NaN,NaN,'l2-sym',1))
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

[d,n]=size(X);
DL2sym=NaN(1,n);
Kn=ones(n,n);
Pn=ones(1,n);
for k=1:d
    Xk=X(k,:);
    if n>1
        Dk=squareform(pdist(Xk'));
    else
        Dk=0;
    end
    Kn=2*Kn.*(-Dk+1);
    pk=2*Xk-2*Xk.^2+1;
    Pn=Pn.*pk;
end
for i=1:n
    DL2sym(i)=sum(sum(Kn(1:i,1:i)))/i^2-2*sum(Pn(1:i))/i;
end
DL2sym=DL2sym+(4/3)^d;
DL2sym=sqrt(DL2sym);
end

    
