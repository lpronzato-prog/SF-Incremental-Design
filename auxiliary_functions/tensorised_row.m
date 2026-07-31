function [ MT ] = tensorised_row( row )
% function [ MT ] = tensorised_row( row,d )
%__________________________________________________________________________  
% row = [k_1,k_2,...,k_m] is a row vector of length m, with sum_i k_i=d
% For example, row = [2    0     1     1     0    0 ] with d=4
% MT is a mm*d matrix, with mm = multinomial(k_1,...,k_m)=d!/(k_1!...k_m!)
% indicating all possible rows that can be constructed from values a_1,...a_m, 
% when a_i can be chosen k_i times
%-----
% L. Pronzato, 2025
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
%-----

m=length(row);
d=sum(row);
I_nonzero=find(row~=0);
LI_nonzero=length(I_nonzero);

if LI_nonzero==1
    MT=I_nonzero(1)*ones(1,d);
else    
    MT=[];
    for i=1:LI_nonzero
        ii=I_nonzero(i);
        rowi=row; rowi(ii)=rowi(ii)-1;
        LTi=log(factorial(d-1))-sum(log(factorial(rowi))); Ti=round(exp(LTi));
    
        Tens_rowi=tensorised_row( rowi );
        MT=[MT
            ii*ones(Ti,1) Tens_rowi];
    end
end

end

