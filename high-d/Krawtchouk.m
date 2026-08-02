function [ P ] = Krawtchouk( j,x,n,s )
% function [ P ] = Krawtchouk( j,x,n,s )
%__________________________________________________________________________
% computes the value of the Krawtchouk polynomial P_j(x;n,s)
% P_j(x;n,s) = sum_i=0^j (-1)^i*(s-1)^(j-i)*nchoosek(x,i)*nchoosek(n-x,j-i)
% All variables must be scalar
%-----
% Author: L. Pronzato, 2020 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE.

P=0;
for i=0:1:j
    P=P+(-1)^i*(s-1)^(j-i)*(gamma(x+1)/gamma(i+1)/gamma(x-i+1))*...
        (gamma(n-x+1)/gamma(j-i+1)/gamma(n-x-j+i+1));
end

end
