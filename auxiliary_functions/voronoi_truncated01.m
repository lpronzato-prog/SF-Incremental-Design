function [ dmax, Dmax, Candidates, i_Candidates ] = voronoi_truncated01( X, A, b, show )
% function [ dmax, Dmax, Candidates, i_Candidates ] = voronoi_truncated01( X, A, b, show )
%
% X = design (d = dimension * n = nb. of points) in [0,1]^d
% A, b define m additional linear constraints Ax<= b (A is m*d, b is m*1)
% to the hypercube [0,1]^d (all point in X should be feasible!)
% show==1 to plot Voronoi diagram when d=2
% Uses polytop_new.m
%
% dmax = max distance between a point of [0,1]^d and a point of X
% Dmax = 1*n vector of max distances to points of X
% Candidates = matrices d * nc of coordinates of points at max distance
%   from points in X in Voronoi cells
% i_Candidates(k) = i indicates that points Candidates(:,k) are at max
%   distance from X(:,i) in its Voronoi cell
% 
% Use [ dmax, ~ ] = voronoi_truncated01( X, [], [], 0 ) to commpute only minimax
% distance criterion in [0,1]^d
%-----
% Author: L. Pronzato, 2015 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

Xp=unique(X','rows');
X=Xp';

[d,~]=size(X);
factor=sqrt(d)+1/2;
toler=1e-6;

% 1) Create XX=X + some points far away so that the voronoi diagram close
% to [0,1]^d is not modified
XX=[X ones(d,d)/2+factor*eye(d) ones(d,d)/2-factor*eye(d)];

if d==2 && show==1
    plot(XX(1,:),XX(2,:),'rp','MarkerSize',5)
end    

[V,C] = voronoin(XX');
nc=0;
i_Candidates=[];
Candidates=[];
Dmax=[];
nv_tot=0;
for i=1:length(C)
     if all(C{i}~=1)
        % No vertex at infinity
        % X(:,i) is center of Voronoi cell
        Vertices = V(C{i},:);
        Cell = convhulln(Vertices);
        [nve,~]=size(Vertices);
        [nfaces,d]=size(Cell); 
        H=zeros(nve,nfaces);
        for ii=1:nfaces
            H(Cell(ii,:)',ii)=1;
        end
        for k=1:d
            ek=zeros(1,d); ek(k)=1; 
            nek=zeros(1,d); nek(k)=-1;
            [Vertices,H]=polytop_new(ek,1,Vertices,H,toler);
            [Vertices,H]=polytop_new(nek,0,Vertices,H,toler);
        end
        if isempty(b)==0
            % There are additional constraints Ax<=b
            for ic=1:length(b)
                [Vertices,H]=polytop_new(A(ic,:),b(ic),Vertices,H,toler);
            end
        end

        if d==2 && show==1
            CC=convhull(Vertices(:,1),Vertices(:,2));
            patch(Vertices(CC,1),Vertices(CC,2),i)
            voronoi(X(1,:),X(2,:),'k--')
        end
       
        % find vertices at max distance from center X(:,i)
        Vertices=Vertices';
        [~,nv]=size(Vertices);
        nv_tot=nv_tot+nv;
        Norms_2=sum((Vertices-X(:,i)*ones(1,nv)).^2);
        dmax_2=max(Norms_2); 
        imax=find(Norms_2==dmax_2);
        Candidates=[Candidates Vertices(:,imax)];
        i_Candidates=[i_Candidates i*ones(1,length(imax))];
        Dmax=[Dmax sqrt(dmax_2)];
        nc=nc+length(imax);
     end
end
dmax=max(Dmax);

end

