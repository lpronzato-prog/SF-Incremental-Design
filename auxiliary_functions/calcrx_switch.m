function rx=calcrx_switch(X,theta,p,x,cov_type,TENSORISED)
% function rx=calcrx_switch(X,theta,p,x,cov_type,TENSORISED) 
%-------------------
% For K a kernel defined by cov_type (e.g., 'Matern32', or 'Riesz', or 
%   'Gaussian'), with correlation length 1/theta (with p an additional
%   parameter in K, such as the power for a Riesz kernel), 
% X (d*n) an n-point design in R^d, 
% x (d*m) an m-point design in R^d,
% computes the corresponding n*m kernel matrix rx
% K is isotropic, unless TENSORISED = 1 for which K is a tensor product
%-----
% L. Pronzato, 2025
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
%-----

[d,n]=size(X);
[~,m]=size(x);

if TENSORISED==1 
    % tensorised kernel
    % compute the product of matrices, component by component
    rx=ones(n,m);
    for i=1:d
        rxi=calcrx_switch(X(i,:),theta,p,x(i,:),cov_type,0);
        rx=rx.*rxi;
    end
else % not a tensorised kernel 
    D=pdist2(X',x');

    switch lower(cov_type)
        case 'gaussian'
            rx=exp(-theta*D.^2);
        case 'matern12'
            rx= exp(-theta*D);    
        case 'matern32'
            rx= (1+sqrt(3)*theta*D).*exp(-sqrt(3)*theta*D);
        case 'matern32p1' % (1+K3/2)
            rx= (1+sqrt(3)*theta*D).*exp(-sqrt(3)*theta*D)+1;                        
        case 'matern52'
            rx=(1+sqrt(5)*theta*D+theta^2*D.^2*5/3).*exp(-sqrt(5)*theta*D);
        case 'distance-short'  % distance metric 1-theta*||x-x'||^p
            rx=-theta*D.^p +1;
        case 'multiquadric'  % ~ (||x-x'||^2*theta^2+1)^(-p)
            rx=(theta^2*D.^2+1).^(-p);
        case 'riesz'  % Riesz kernel (||x-x'||^(-p), 0<p<1)
            rx=D.^(-p);
        case 'triangular'  % triangular kernel
            rx=max(-theta*D +1,0);  
        case 'cosine'  % cosine kernel
            rx=(cos(theta*pi*D));                       
        case 'l2-wa' % kernel of discrepancy_wa
            rx=-D+D.^2+3/2;
        case 'l2-sym'  % kernel of discrepancy_L2sym
            rx=2*(-D+1);
        case 'l2-cent'  % kernel of discrepancy_cent
            rx=-D/2+1+abs(X'-1/2)*ones(1,m)/2+ones(n,1)*abs(x-1/2)/2;
        case 'l2-star'  % kernel of discrepancy_L2star
            Dmax=max(X'*ones(1,m),ones(n,1)*x);
            rx=-Dmax+1;
        case 'l2-starm'  % kernel of discrepancy_L2starM
            Dmax=max(X'*ones(1,m),ones(n,1)*x);
            rx=-Dmax+2;
        case 'l2-extreme'  % kernel of discrepancy_L2extreme
            Dmax=max(X'*ones(1,m),ones(n,1)*x);
            Dmin=min(X'*ones(1,m),ones(n,1)*x);
            rx=(-Dmax+1).*Dmin;
        case 'sinc'  % sinc kernel 
            rx=sin(theta*D)./(theta*D);
    end    
end    
end
