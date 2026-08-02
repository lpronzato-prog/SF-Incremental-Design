function [aN,bN]=GoldenS(f,a0,b0,epsi)
% function [aN,bN]=GoldenS(f,a0,b0,delta)
%__________________________________________________________________________
% Golden Section method for minimization of f on [a0,b0]
% returns [aN,bN], with bN-aN<epsi, which contains a minimizer
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

a=a0; b=b0;
delta=b-a; 
if delta<epsi
    aN=a0; bN=b0;
else    
    phi=(sqrt(5)-1)/2;  % GS
    e=a+phi*delta; fe=feval(f,e);
    while delta>=epsi
        if e < (a+b)/2
            ep=a+phi*delta; 
            fu=fe; fv=feval(f,ep); u=e; v=ep;
        else
            ep=a+(1-phi)*delta; 
            fu=feval(f,ep); fv=fe; u=ep; v=e;
        end
        if fu<=fv
            b=v; e=u; fe=fu;
        else
            a=u; e=v; fe=fv;
        end
        delta=phi*delta;
    end
    aN=a; bN=b;
end    

