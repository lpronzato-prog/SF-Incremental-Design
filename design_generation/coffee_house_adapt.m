function [ X,Time ] = coffee_house_adapt( Xcand, nmax, Alpha, norepeat, center)
% function [ X,Time ] = coffee_house_adapt( Xcand, nmax, Alpha, norepeat, center)
%__________________________________________________________________________
% = Algorithm 6.1: Greedy packing 
%       of [Karvonen, Pronzato and Zhigljavsky, 2026]
% The design space is [0,1]^d. If Alpha=zeros(1,nmax) (no boundary 
% avoidance) and center=0 (initilization at Xcand(1)), the algorithm  can be 
% applied to any design space XX, provided that Xcand contains points in XX.  
%__________________________________________________________________________
% = greedy packing (with boundary avoidance when Alpha(n)>0, n=1,...,nmax) 
% computes a design X (d*nmax) by boundary-phobic greedy packing 
% (The indices of its points in Xcand are passed as global variable IndX_g
%   if needed)
% The design X is searched within Xcand (d*Q), a set of Q candidate points
% No computation of Q*Q interdistances matrix is involved
% alpha = Alpha(n) at iteration n defines the amount of bounday avoidance:
%   it corresponds to 1/beta in [Nogalez-Gómes, Pronzato and Rendas, 2021], 
%   where beta=d/(2*(nmax*Vd)^(-1/d))-sqrt(d). 
%   The value beta=2*sqrt(2*d) is recommended in [Shang and Apley 2020].
% Alpha = zeros(1,nmax) gives the standard greedy-packing 
%   (= coffee-house algorithm)
%----- input variables
%   Xcand (d*Q) = a set of Q >= nmax candidate points in [0,1]^d 
%   nmax = number of points generated 
%   Alpha (1*nmax) = relaxations coefficients that set the level of bounday
%                       avoidance
%   norepeat = 1 to force no repetitions of design points (which may happen when alpha>0, 
%                   e.g., if all remaining candidates are on the boundary)
%   center = 1 to initialize x(1) at the point in Xcand closest to the
%                   center of the cube [0,1]^d, at Xcand(1) otherwise
%----- output variables
% X (d*nmax) = the design generated
% Time (1*nmax) = running time along iterations    
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

[d,Q]=size(Xcand);
t1=tic;
Time=NaN(1,nmax);  

IndX=zeros(1,nmax);

%------------------
% 1) Initialization
%------------------
Xall=NaN(d,nmax);
if center==1
    % initialize by "central point" (in case of a design [0,1]^d) 
    Center=ones(d,1)/2; 
    D_to_center=pdist2(Center',Xcand');
    [~,i1]=min(D_to_center); 
else 
    i1=1;
end    

X=Xcand(:,i1); IndX(1)=i1;
Xall(:,1)=X;
D_to_X=pdist2(Xcand',X');  % distances to design points

alpha1=Alpha(1);
Xboundary=[Xcand' -Xcand'+1]/alpha1;

%--------------------
% 2) Start iterations
%--------------------
timebar= waitbar(0,'Coffee-house...'); 
Time(1)=toc(t1);
for n=2:nmax
    waitbar(n/nmax,timebar);      
    alpha=Alpha(n);
    if alpha~=alpha1 % useless to repeat always the same calculation otherwise
        Xboundary=[Xcand' -Xcand'+1]/alpha;
    end
    [distmax,~]=max(D_to_X);  
    Imax=find(D_to_X==distmax);
    %Imax=find(abs(D_to_X-distmax)<1e-1);
    
    ibest=Imax(1);
    if alpha>0
        if norepeat==0 
            [D_to_X Xboundary]
            [~,ibest]=max(min([D_to_X Xboundary],[],2));
            ibest=ibest(1);
        else
            Ind_free=setdiff(1:Q,IndX);
            Dum=[D_to_X Xboundary];
            Dum=Dum(Ind_free,:);
            [~,ibest]=max(min(Dum,[],2)); ibest=ibest(1);
            ibest=Ind_free(ibest);
        end
    end
    Xnew=Xcand(:,ibest);
    D_to_X=min(D_to_X,pdist2(Xcand',Xnew')); % update all D_to_X     
    x_last=Xnew; 
    Xall(:,n)=x_last;
    X=Xall(:,1:n);
    IndX(n)=ibest;  
    
    Time(n)=toc(t1);
end
close(timebar) 
end
        