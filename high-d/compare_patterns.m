function [ ibest ] = compare_patterns( B1,B2 )
% function [ ibest ] = compare_patterns( B1,B2 )
%__________________________________________________________________________
% Compares two distance patterns (or word-lenght patters) B1 and B2 (row
% vectors of length d+1)
% Returns 1 if B1 is better than B2 (that is B1(k)<B2(k) for k the fist index
% where they differ) and 2 otherwise
%-----
% Author: L. Pronzato, 2020 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE.

Diff=B1-B2;
k=find((Diff~=0),1);
ibest=1;
if B2(k)<B1(k), ibest=2; end

end

