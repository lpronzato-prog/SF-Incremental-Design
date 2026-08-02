function Deltaopt = Deltaopt_cube(alpha,d,n,s,accuracy)
% function Deltaopt = Deltaopt_cube(alpha,d,n,s,accuracy)
%__________________________________________________________________________
% Computes the optimal size of a cube [-Delta,Delta]^d that minimizes the 
% (approximated) expected distortion ( = (Ls-mean quantization error)^s) 
% for a random design with n points independently distributed with the 
% spherically symmetric beta distribution in [-Delta,Delta]^d computed 
% by density_norm_X(rho,alpha,Delta,d),
% for the reference measure mu uniform in the cube [-1,1]^d.
% accuracy = precision required on Ropt (e.g., 1e-4)
% Corresponds to Section 7.2 of the book [Karvonen, Pronzato and Zhigljavsky, 
% Springer, 2026] 
% When alpha = 1, the design points are uniform in [-Delta,Delta]^d
% when alpha = 0 and n<2^d, the design points are randomly sampled without
%   replacement among the vertices of the cube [-Delta,Delta]^d. One may
%   altenatively sample the vertices of a fractional-factorial selection of
%   these vertices to ensure a large packing radius (see 
%   Example_fractional_factorial.m for an example), or apply the
%   greedy-packing algorithm using this fractional-factorial selection of 
%   vertices as candidate set (see Examples_generation_cube.m).
%
% Uses Edistortion_cube.m
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

% E{quantization error^s} for n points uniformly distributed in the ball or on the sphere

ED_Delta=@(delta) Edistortion_cube(alpha,delta,d,n,s);
[aN,bN]=GoldenS(ED_Delta,0.1,1,accuracy);  
Deltaopt=(aN+bN)/2;
end
