function [ PotK2 ] = potential_K2( x, theta,p,cov_type,TENSORISED)
% function [ PotK2 ] = potential_K2( x, theta,p,cov_type,TENSORISED)
%-------------------
% For dimension d>1, REQUIRES TENSORISED = 1,
%   only for cov_type = 'Matern12', 'Matern32' or 'Gaussian'  
%-------------------
% For K a kernel defined by cov_type (e.g., 'Matern32'), 
% with correlation length 1/theta (the additional parameter p is not used),
% computes PotK2(x) = int_0^1 K^2(x,t) dt 
% x (d*n) must have its components in [0,1] 
% returns PotK2 which is 1*n
%-----
% L. Pronzato, 2025
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
%-----

% computes s(x,y)=int_0^1 K^2(x,t) dt for a tensorised kernel (i.e., TENSORISED=1)
% Matern12, Matern32, Gaussian

[d,n]=size(x);
PotK2=ones(1,n);

if d>1 && TENSORISED==0
    PotK2=zeros(1,n);
    display('Error: K must be a tensorised kernel')
    return
end
if TENSORISED==1 
    for i=1:d
        [ PotK2i ] = potential_K2( x(i,:), theta,p,cov_type,0);
        PotK2=PotK2.*PotK2i;
    end
else
    switch lower(cov_type)
        case 'gaussian' % = potential for K but 2*theta
            [ PotK2, ~ ] = potential_energy( x, 2*theta,p,cov_type,0);
        case 'matern12' % = potential for K but 2*theta
            [ PotK2, ~ ] = potential_energy( x, 2*theta,p,cov_type,0);   
        case 'matern32' % need to compute...
            T=theta*sqrt(3);
            ineg=find(x<0);
            igt1=find(x>1);
            iok=find((x>=0) & (x<=1));
            if isempty(ineg)==0
                z=x(ineg);
                PotK2(ineg)= -(-5*exp(2*T*z)+6*exp(2*T*z).*z*T-2*T^2 ...
                *exp(2*T*z).*z.^2+5*exp(2*(z-1)*T)-6*T*exp(2*(z-1)*T).*z...
                    +2*T^2*exp(2*(z-1)*T).*z.^2-4*T^2*exp(2*(z-1)*T).*z...
                    +6*T*exp(2*(z-1)*T)+2*T^2*exp(2*(z-1)*T))/T/4;
            end
            if isempty(iok)==0
                z=x(iok);
                PotK2(iok)= -(5*exp(2*T)+6*exp(2*T)*T*z+2*exp(2*T)*...
                    T^2*z.^2+5*exp(4*T*z)+6*exp(4*T*z)*T-6*exp(4*T*z)*...
                    T.*z-10*exp(2*T*z+2*T)+2*exp(4*T*z)*T^2.*z.^2-...
                    4*exp(4*T*z)*T^2.*z+2*exp(4*T*z)*T^2).*...
                    exp(-2*T*z-2*T)/T/4;
            end
            if isempty(igt1)==0
                z=x(igt1);
                PotK2(igt1)=-(5*exp(-2*T*z)+6*exp(-2*T*z)*T.*z+2*...
                    exp(-2*T*z)*T^2.*z.^2-5*exp(-2*(z-1)*T)-6*...
                    exp(-2*(z-1)*T)*T.*z-2*exp(-2*(z-1)*T)*T^2.*z.^2+...
                    4*exp(-2*(z-1)*T)*T^2.*z+6*exp(-2*(z-1)*T)*T-2*...
                    exp(-2*(z-1)*T)*T^2)/T/4;
            end
    end
end    

end

