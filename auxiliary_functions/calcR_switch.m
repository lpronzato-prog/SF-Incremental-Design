function R=calcR_switch(X,theta,p,cov_type,TENSORISED)
% function R=calcR_switch(X,theta,p,cov_type,TENSORISED) 
%-------------------
% For K a kernel defined by cov_type (e.g., 'Matern32', or 'Riesz', or 
%   'Gaussian'), with correlation length 1/theta (with p an additional
%   parameter in K, such as the power for a Riesz kernel), computes the
%   kernel matrix R for the n design points in X (a d*n matrix).
% K is isotropic, unless TENSORISED = 1 for which K is a tensor product
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 
%-----

% n*n distance matrix
[d,n]=size(X);

if TENSORISED==1 
    % tensorised kernel
    % compute the product of matrices, component by component
    R=ones(n,n);
    for i=1:d
        Ri=calcR_switch(X(i,:),theta,p,cov_type,0);
        R=R.*Ri;
    end
else % not a tensorised kernel   
    if n>1
        D=squareform(pdist(X'));
    else
        D=0;
    end

    switch lower(cov_type)
        case 'gaussian'
            R=exp(-theta*D.^2);
        case 'matern12'
            R= exp(-theta*D);
        case 'matern32'
            R= (1+sqrt(3)*theta*D).*exp(-sqrt(3)*theta*D);
        case 'matern32p1' % (1+K3/2)
            R= (1+sqrt(3)*theta*D).*exp(-sqrt(3)*theta*D) +1;                       
        case 'matern52'
            R=(1+sqrt(5)*theta*D+theta^2*D.^2*5/3).*exp(-sqrt(5)*theta*D);
        case 'distance-short'  % distance metric 1-theta*||x-x'||^p
            R=-theta*D.^p +1;
        case 'multiquadric'  % ~ (||x-x'||^2*theta^2+1)^(-p)
            R=(theta^2*D.^2+1).^(-p); 
        case 'riesz'  % Riesz kernel (||x-x'||^(-p), 0<p<1)
            R=D.^(-p);
        case 'triangular'  % triangular kernel
            R=max(-theta*D +1,0);
        case 'cosine'  % cosine kernel
            R=(cos(theta*pi*D));            
        case 'l2-wa' % kernel of discrepancy_wa
            R=-D+D.^2+3/2;
        case 'l2-sym'  % kernel of discrepancy_L2sym
            R=2*(-D+1);
        case 'l2-cent'  % kernel of discrepancy_cent
            R=-D/2+1+abs(X'-1/2)*ones(1,n)/2+ones(n,1)*abs(X-1/2)/2;
        case 'l2-star'  % kernel of discrepancy_L2star
            Dmax=max(X'*ones(1,n),ones(n,1)*X);
            R=-Dmax+1;
        case 'l2-starm'  % kernel of discrepancy_L2starM
            Dmax=max(X'*ones(1,n),ones(n,1)*X);
            R=-Dmax+2;
        case 'l2-extreme'  % kernel of discrepancy_L2extreme
            Dmax=max(X'*ones(1,n),ones(n,1)*X);
            Dmin=min(X'*ones(1,n),ones(n,1)*X);
            R=(-Dmax+1).*Dmin;          
        case 'sinc'  % sinc kernel 
            R=sin(theta*D)./(theta*D);
            R=triu(R,1)+tril(R,-1)+eye(n); % to replace NaN on the diagonal by 1
    end
end
end
