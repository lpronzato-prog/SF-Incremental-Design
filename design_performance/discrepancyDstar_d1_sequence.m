function [ discr_star_n, discr_n ] = discrepancyDstar_d1_sequence( x )
% function [ discr ] = discrepancyDstar_d1_sequence( x )
%__________________________________________________________________________
% computes the star and extreme discrepancies of the sequence
% x(1),...,x(n); 
% discr_star_n, discr_n are vectors of length n
%-----
% Author: L. Pronzato, 2022 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE.

n=length(x);
discr_star_n=NaN(1,n);
discr_n=NaN(1,n);
for i=1:n
    [discr_star_n(i), discr_n(i)]=discrepancyDstar_d1( x(1:i) );
end

