function [ Xcand ] = candidate_set( d,q,k,center,Scramble )
% function [ Xcand ] = candidate_set( d,q,k,center,Scramble )
%-------------------
% Generates a candidate set Xcand 
% If q>0, Xcand = Sobol'(d,2^q) complemented by a k^d
%   fractional factorial design if k>=2
% If q=0, Xcand is a k^d fractional factorial design
% If center=1, Xcand contains the center 0.5*ones(d,1) of the cube
% If Scramble=1, the Sobol' set is scrambled by scramble(pS,'MatousekAffineOwen');
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

if q==0
    Xcand=[];
else    
    pS = sobolset(d);
    if Scramble==1
        pS = scramble(pS,'MatousekAffineOwen');
    end   
    Q=2^q;
    Xcand=(pS(1:Q,:))';
end
if k>=2
    Complement=((fullfact(k*ones(1,d))-1)/(k-1))'; % for any d
    %size(Complement)
    Xcand=[Xcand Complement]; 
end    
if center==1
    Xcand=[1/2*ones(d,1) Xcand ]; 
end
Xcand=unique(Xcand','rows','stable'); Xcand=Xcand';

end

