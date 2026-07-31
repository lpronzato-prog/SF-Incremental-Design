function [s,H]=polytop_new(XX,c,s,H,toler)
% function [s,H]=polytop_new(XX,c,s,H,toler)
% -------------------------------------------------------------------------
% Computes the intersection of a polytope P defined by its vertices s with 
% the half-space XX*theta <= c 
% XX = row vector (1*p)
% c = scalar
% s = a matrix (m*p) whose each row contains the coordinates of a vertex of P
% H = a matrix that summarizes the adjacency relationships between the faces and vertices of P
% toler = tolerance for vertex merging (new vertex s=a*s1+(1-a)*s2,
% with s1 (resp. s2) an old exterior (resp. intérior) vertex, s is confused
% with s2 when a<toler) and for vertex construction (s exterior if XX*s>c+toler).
%-----
% Author: L. Pronzato, 1992 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

p=length(XX); 
	[ns,p]=size(s);
	% nø des sommets à l'extérieur de la contrainte suivante
	is=find(XX(1,:)*s'>c+toler);
	ls=length(is);
	if ls==ns
		% tous les sommets sont à l'extérieur
		text='polyèdre vide'
		s=[];
		return
	elseif ls~=0
		% il existe des sommets à l'extérieur
		% partition des sommets
		He=H(is,:); se=s(is,:);	% sommets extérieurs
		H(is,:)=[]; s(is,:)=[];	% sommets intérieurs
		csH=zeros(ns-ls,1);	% colonne supplémentaire pour H
		% construction des nouveaux sommets
		snew=[]; Hnew=[];
		for j=1:ls
			% définition des sommets intérieurs adjacents au
			% jème sommet extérieur
			iadj=find(H*He(j,:)'>=p-1);
			ladj=length(iadj);
			if ladj~=0
			    for k=1:ladj
                    % test supplémentaire d'adjacence
                    kadj=iadj(k);
                    HHe=H(kadj,:).*He(j,:);
                    kki=find(H*HHe'>=p-1);
                    kke=find(He*HHe'>=p-1);
                    if length(kki)==1 & length(kke)==1
                        % pas d'autre sommet sur l'arête: il y
                        % a bien adjacence
                        lam=(XX*s(kadj,:)'-c)/(XX*(s(kadj,:)-se(j,:))');
                        if lam<=toler % nouveau sommet=ancien 
				    		csH(kadj)=1; % ancien sommet sur la nouvelle contrainte i
                        else
				       	    snew=[snew;(1-lam)*s(kadj,:)+lam*se(j,:)];
				    		Hnew=[Hnew;He(j,:).*H(kadj,:) 1];
				    		% une colonne de plus pour la contrainte i
                        end
                    end
			   end
			end
		end
		% mise à jour du tableau H
		% une colonne supplémentaire au tableau des sommets intérieurs
		H=[H csH];
		% fusion des tableaux
		H=[H;Hnew]; 
		s=[s;snew];
		% élimination des hyperplans inutiles (contenant moins de p
		% sommets)
		v=(sum(H)<p);
		H(:,v)=[]; 
	end
