function [ discr_star_n, discr_n ] = discrepancyDstar_d1( x )
% function [ discr ] = discrepancyDstar_d1( x )
%__________________________________________________________________________
% computes the star and extreme discrepancies of x (vector of length n)
%-----
% Author: L. Pronzato, 2022 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE.

n=length(x);
x=x(:); x=x';
[xsort,~]=sort(x);
discr_star_n = 1/(2*n)+max(abs(xsort-(2*(1:n)-1)/(2*n)));

Deltan=(1:n)/n-xsort;
discr_n = 1/n+max(Deltan)-min(Deltan);
end

