function [ PR ] = Packing_Radius( Xn, sequence )
% function [ PR ] = Packing_Radius( Xn, sequence )
%__________________________________________________________________________
% For Xn (d*n) an n-point design, returns 
% PR of size m = length(sequence), a row vector containing the packing 
% radii of the sequence of designs Xn(:,(1:i1)),...,Xn(:,(1:im)), with 
% [i_1 ... i_m] = sequence = a vector of m consecutive integers in {1,...,n}
% If sequence = k, PR is just PR(Xn(:,(1:k)))
%-----
% Author: L. Pronzato, 2021 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

[d,n]=size(Xn);

LPR=length(sequence);
PR=NaN(1,LPR);
for i=1:LPR
    iCR=sequence(i);
    if iCR>1 % no packing radius for a one point design
        PR(i)=min(pdist(Xn(:,1:iCR)'));
    end 
end
PR=PR/2;
end
