function [ WLP ] = From_distancepattern_to_WLP( Dist_pattern,d,m )
% function [ WLP ] = From_distancepattern_to_WLP( Dist_pattern,d,m )
%__________________________________________________________________________
%   Computes the Word Length Pattern WLP for a 2^d-m fractional design with
%   given distance pattern
%-----
% Author: L. Pronzato, 2020 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE.

k=d-m;
WLP=zeros(1,d+1);
for jj=1:d+1   % +1 since we need to shift indices (stupid matlab)  
    WLPj=0;
    for ii=1:d+1; % +1 since we need to shift indices (stupid matlab)
        i=ii-1;
        j=jj-1;
        WLPj=WLPj+Dist_pattern(ii)*Krawtchouk( j,i,d,2 );
    end
    WLP(jj)=WLPj/2^k;
end


end

