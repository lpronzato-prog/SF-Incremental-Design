function M=tri2p01(p)
% function M=tri2p01(p)
% -----------------
% For list1 and list2 two lists of p objects (i.e., vectors of size p),
% returns M, a (2^p)*2p matrix of 0 and 1, indicating all lists of p
% objects that can be taken from list1 or list2. 
% Example: the kth row [v1 v2] of M is such that v1(i)=1 and v2(i)=0, for 
% i=1,...,p, if the ith object of the kth list is taken from list1.
%--------------
% called by coffee_house_subspaces.m
%-----
% L. Pronzato, 2025
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
%-----

if p==1
	M=eye(2);
else
	MM=tri2p01(p-1);
	M=[ones(2^(p-1),1) MM(:,1:p-1) zeros(2^(p-1),1) MM(:,p:2*p-2);...
	   zeros(2^(p-1),1) MM(:,1:p-1) ones(2^(p-1),1) MM(:,p:2*p-2)];
end
