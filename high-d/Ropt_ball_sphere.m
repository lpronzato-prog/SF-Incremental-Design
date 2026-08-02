function Ropt = Ropt_ball_sphere(set,d,n,s,accuracy)
% function Ropt = Ropt_ball_sphere(set,d,n,s,accuracy)
%__________________________________________________________________________
% Computes the optimal radius Ropt that minimizes the expected distortion 
% ( = (Ls-mean quantization error)^s) for a random design with n points 
% independently uniformly distributed 
%   - in the ball B_d(0,R) when set = 'ball'
%   - on the sphere S_{d-1}(0,R) when set = 'sphere'
% for the reference measure mu uniform in B_d(0,1).
% accuracy = precision required on Ropt (e.g., 1e-4)
% Corresponds to Section 7.1.4 of the book [Karvonen, Pronzato and Zhigljavsky, 
% Springer, 2026] 
%
% Uses Edistortion_ball.m, Edistortion_sphere.m
%-----
% Author: L. Pronzato, 2025 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE. 

% E{quantization error^s} for n points uniformly distributed in the ball or on the sphere
switch lower(set)
    case 'ball'
        ED_R=@(R) Edistortion_ball(R,d,n,s);
    case 'sphere'
        ED_R=@(R) Edistortion_sphere(R,d,n,s);
    end    
[aN,bN]=GoldenS(ED_R,0.1,1,accuracy);  
Ropt=(aN+bN)/2;
end
