function [ discr_star_n, discr_n ] = star_extreme_radial_discr_sequence( X )
% function [ discr_star_n, discr_n ] = star_extreme_radial_discr_sequence( X )
%__________________________________________________________________________
% For X (d*n) an n-point design in [0,1]^d, returns discr_star_n, discr_n, 
% the 1*n vectors of the star and extreme discrepancies of radial distances
% max_k |X(k,i)-1/2|, for i=1,...,n
% Uses the property that for U uniformly distributed in [0,1]^d, 
% z=(2*max_i |U_i-1/2|)^d is uniformly distributed in [0,1]
% Uses discrepancyDstar_d1_sequence.m
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

[d,~]=size(X);
D_radial=max(abs(X-1/2),[],1);
[ discr_star_n, discr_n ] = discrepancyDstar_d1_sequence( (2*D_radial).^d );

end

