function [ kbarx, Kbar ] = potential_energy( x, theta,p,kernel,TENSORISED)
% function [ kbarx, Kbar ] = potential_energy( x, theta,p,kernel,TENSORISED)
%-------------------
% For dimension d>1, REQUIRES TENSORISED = 1,
%   not all kernels of calcR_switch.m can be used
%-------------------
% For K a kernel defined by cov_type (e.g., 'Matern32', or 'Riesz', or 
%   'Gaussian'), with correlation length 1/theta (with p an additional
%   parameter in K, such as the power for a Riesz kernel), 
% computes kbarx = int_[0,1]^d K(x,t) dt and Kbar = int_[0,1]^d int_[0,1]^d K(x,t) dt dx
%-----
% L. Pronzato, 2025
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
%-----

[d,n]=size(x);

if d>1 && TENSORISED ~= 1 
    disp('Only works for tensorised kernels')
end    
if d>1
    % use a tensorised kernel
    % compute the product of matrices, component by component
    Kbar=1;
    kbarx=ones(1,n);
    for i=1:d
        [ kbarxi, Kbari ] = potential_energy( x(i,:), theta,p,kernel,0);
        Kbar=Kbar*Kbari;
        kbarx=kbarx.*kbarxi;
    end
else % not a tensorised kernel 

    switch lower(kernel)
        case 'gaussian'
            theta=abs(theta);
            kbarx=(1/2)*sqrt(pi)*(erf(sqrt(theta)*x)-erf(-sqrt(theta)+sqrt(theta)*x))/sqrt(theta);
            Kbar=-(exp(theta)*sqrt(theta)-sqrt(theta)-erf(sqrt(theta))*sqrt(pi)*theta*exp(theta))*exp(-theta)/theta^(3/2);
        case 'matern32'
            % Adapted from [Ginsbourger, Roustant, Schuhmacher, Lenz, 2014]
            A1=@(u) (2+u).*exp(-u);
            T=theta*sqrt(3);
            xi1=1/T;
            Kbar=2*xi1*(2-3*xi1+(1+3*xi1)*exp(-1/xi1));
            ineg=find(x<0);
            igt1=find(x>1);
            iok=find((x>=0) & (x<=1));
            if isempty(iok)==0
                kbarx(iok)=xi1*(4-A1(x(iok)/xi1)-A1((1-x(iok))/xi1)); 
            end
            if isempty(ineg)==0
                y=x(ineg);
                kbarx(ineg)=(2*exp(T*y)-exp(T*y)*T.*y-2*exp(T*(y-1))+...
                    exp(T*(y-1))*T.*y-exp(T*(y-1))*T)/T;
            end
            if isempty(igt1)==0
                y=x(igt1);
                kbarx(igt1)=-(2*exp(-T*y)+exp(-T*y)*T.*y-2*exp(-T*(y-1))-...
                    exp(-T*(y-1))*T.*y+exp(-T*(y-1))*T)/T;
            end
        case 'matern32p1' % (1+K3/2)
            % Adapted from [Ginsbourger, Roustant, Schuhmacher, Lenz, 2014]
            A1=@(u) (2+u).*exp(-u);
            T=theta*sqrt(3);
            xi1=1/T;
            Kbar=2*xi1*(2-3*xi1+(1+3*xi1)*exp(-1/xi1));
            ineg=find(x<0);
            igt1=find(x>1);
            iok=find((x>=0) & (x<=1));
            if isempty(iok)==0
                kbarx(iok)=xi1*(4-A1(x(iok)/xi1)-A1((1-x(iok))/xi1)); 
            end
            if isempty(ineg)==0
                y=x(ineg);
                kbarx(ineg)=(2*exp(T*y)-exp(T*y)*T.*y-2*exp(T*(y-1))+...
                    exp(T*(y-1))*T.*y-exp(T*(y-1))*T)/T;
            end
            if isempty(igt1)==0
                y=x(igt1);
                kbarx(igt1)=-(2*exp(-T*y)+exp(-T*y)*T.*y-2*exp(-T*(y-1))-...
                    exp(-T*(y-1))*T.*y+exp(-T*(y-1))*T)/T;
            end
            Kbar=Kbar+1;
            kbarx=kbarx+1;
        case 'matern12'
            Kbar=(2*(-1+exp(-theta)+theta))/theta^2;
            ineg=find(x<0);
            igt1=find(x>1);
            iok=find((x>=0) & (x<=1));
            if isempty(iok)==0
                y=x(iok);
                kbarx(iok)=-(exp(-theta*y)-2+exp(theta*(y-1)))/theta;
            end
            if isempty(ineg)==0
                y=x(ineg);
                kbarx(ineg)=(exp(theta*y)-exp(theta*(y-1)))/theta;
            end
            if isempty(igt1)==0
                y=x(igt1);
                kbarx(igt1)=(-exp(-theta*y)+exp(-theta*(y-1)))/theta;
            end      
        case 'distance-short'  % distance metric 1-theta*||x-x'||^p
            Kbar=-theta/3 +1;
            ineg=find(x<0);
            igt1=find(x>1);
            iok=find((x>=0) & (x<=1));
            if isempty(iok)==0
                y=x(iok);
                kbarx(iok)=1+theta*(-y.^2 + y -1/2);
            end
            if isempty(ineg)==0
                y=x(ineg);
                kbarx(ineg)=1+theta*(y-1/2);
            end
            if isempty(igt1)==0
                y=x(igt1);
                kbarx(igt1)=1+theta*(-y+1/2);
            end
        case 'multiquadric'  % ~ (||x-x'||^2*theta^2+1)^(-1)
            epsi=theta;
            Kbar=(2*epsi*atan(epsi)-log(1+epsi^2))/epsi^2;
            kbarx=(atan(epsi*x)-atan(epsi*(x-1)))/epsi;
        case 'riesz'  % Riesz kernel (||x-x'||^(-p), 0<p<1)
            s=-p;
            Kbar=2/(s^2+3*s+2);
            ineg=find(x<0);
            igt1=find(x>1);
            iok=find((x>=0) & (x<=1));
            if isempty(iok)==0
                y=x(iok);
                kbarx(iok)=((-y+1).^(s+1)+y.^(s+1))/(s+1);
            end
            if isempty(ineg)==0
                y=x(ineg);
                kbarx(ineg)=((-y+1).^(s+1)-(-y).^(s+1))/(s+1);
            end
            if isempty(igt1)==0
                y=x(igt1);
                kbarx(igt1)=(y.^(s+1)-(y-1).^(s+1))/(s+1);
            end
        case 'triangular'  % triangular kernel
            %R=max(-theta*D.^p +1,0);        
            if theta<=1
                Kbar=1-theta/3;
            else
                Kbar=(3*theta-1)/(3*theta^2);
            end
            %
            if 2<=theta
                iA=find(x<-1/theta);
                iB=find((-1/theta<=x) & (x<0));
                iC=find((0<=x) & (x<1/theta));
                iD=find((1/theta<=x) & (x<1-1/theta));
                iE=find((1-1/theta<=x) & (x<1));
                iF=find((1<=x) & (x<1+1/theta));
                iG=find(1+1/theta<=x);
                if isempty(iA)==0
                    y=x(iA);
                    kbarx(iA)=0*y;
                end
                if isempty(iB)==0
                    y=x(iB);
                    kbarx(iB)=theta/2*(y+1/theta).^2;
                end
                if isempty(iC)==0
                    y=x(iC);
                    kbarx(iC)=1/theta-theta/2*(-y+1/theta).^2;
                end
                if isempty(iD)==0
                    y=x(iD);
                    kbarx(iD)=1/theta;
                end
                if isempty(iE)==0
                    y=x(iE);
                    kbarx(iE)=1/theta-theta/2*(y-1+1/theta).^2;
                end
                if isempty(iF)==0
                    y=x(iF);
                    kbarx(iF)=theta/2*(-y+1+1/theta).^2;
                end
                if isempty(iG)==0
                    y=x(iG);
                    kbarx(iG)=0*y;
                end
            elseif 1<=theta && theta<2
                iA=find(x<-1/theta);
                iB=find((-1/theta<=x) & (x<0));
                iC=find((0<=x) & (x<1-1/theta));
                iD=find((1-1/theta<=x) & (x<1/theta));
                iE=find((1/theta<=x) & (x<1));
                iF=find((1<=x) & (x<1+1/theta));
                iG=find(1+1/theta<=x);
                if isempty(iA)==0
                    y=x(iA);
                    kbarx(iA)=0*y;
                end
                if isempty(iB)==0
                    y=x(iB);
                    kbarx(iB)=theta/2*(y+1/theta).^2;
                end
                if isempty(iC)==0
                    y=x(iC);
                    kbarx(iC)=1/theta-theta/2*(-y+1/theta).^2;
                end
                if isempty(iD)==0
                    y=x(iD);
                    kbarx(iD)=1/theta-theta/2*(-y+1/theta).^2-theta/2*(y-1+1/theta).^2;
                end
                if isempty(iE)==0
                    y=x(iE);
                    kbarx(iE)=1/theta-theta/2*(y-1+1/theta).^2;
                end
                if isempty(iF)==0
                    y=x(iF);
                    kbarx(iF)=theta/2*(-y+1+1/theta).^2;
                end
                if isempty(iG)==0
                    y=x(iG);
                    kbarx(iG)=0*y;
                end
            elseif theta<1
                iA=find(x<-1/theta);
                iB=find((-1/theta<=x) & (x<1-1/theta));
                iC=find((1-1/theta<=x) & (x<0));
                iD=find((0<=x) & (x<1));
                iE=find((1<=x) & (x<1/theta));
                iF=find((1/theta<=x) & (x<1+1/theta));
                iG=find(1+1/theta<=x);
                if isempty(iA)==0
                    y=x(iA);
                    kbarx(iA)=0*y;
                end
                if isempty(iB)==0
                    y=x(iB);
                    kbarx(iB)=theta/2*(y+1/theta).^2;
                end
                if isempty(iC)==0
                    y=x(iC);
                    kbarx(iC)=theta*y-theta/2+1;
                end
                if isempty(iD)==0
                    y=x(iD);
                    kbarx(iD)=1+theta*(-y.^2 + y -1/2);
                end
                if isempty(iE)==0
                    y=x(iE);
                    kbarx(iE)=-theta*y+theta/2+1;
                end
                if isempty(iF)==0
                    y=x(iF);
                    kbarx(iF)=theta/2*(-y+1+1/theta).^2;
                end
                if isempty(iG)==0
                    y=x(iG);
                    kbarx(iG)=0*y;
                end            
            end 
        case 'l2-wa' % kernel of discrepancy_wa
            kbarx=4/3;
            Kbar=4/3;
        case 'l2-sym'  % kernel of discrepancy_L2sym
            kbarx=2*x-2*x.^2+1;
            Kbar=4/3;
        case 'l2-cent'  % kernel of discrepancy_cent
            kbarx=0.5*abs(x-1/2)-0.5*(x-1/2).^2+1;
            Kbar=13/12;
        case 'l2-star'  % kernel of discrepancy_L2star
            kbarx=(-x.^2+1)/2;
            Kbar=1/3;
        case 'l2-starm'  % kernel of discrepancy_L2starM
            kbarx=(-x.^2+3)/2;
            Kbar=4/3;
        case 'l2-extreme'  % kernel of discrepancy_L2extreme
            kbarx=x.*(-x+1)/2;
            Kbar=1/12;   
     end 
end
end

