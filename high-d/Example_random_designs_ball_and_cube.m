% Example of construction of random (or quasi-random) designs in high-dimensional balls and cubes
% with small Ls-mean quantization error
% see Chap. 7 of the book [Karvonen, Pronzato and Zhigljavsky, Springer, 2026]

d=30; 
n=1000;
s=2;
accuracy=1e-2;

% The reference measure uniform in the ball B_d(0,1)
    % For n random points uniformly distributed in B_d(0,R), the optimal radius is:
        Ropt = Ropt_ball_sphere('ball',d,n,s,accuracy)
    % the expected distortion (distortion = (Ls-mean quantization error)^s) is:
        Edist = Edistortion_ball(Ropt,d,n,s)
    %
    % For n random points uniformly distributed on S_{d-1}(0,R), the optimal radius is:
        Ropt = Ropt_ball_sphere('sphere',d,n,s,accuracy)
    % the expected distortion (distortion = (Ls-mean quantization error)^s) is:
        Edist = Edistortion_sphere(Ropt,d,n,s)

%% The reference measure uniform in the cube [-1,1]^d
    % For n random points uniformly distributed in [-delta,delta]^d, the optimal delta is:
        Deltaopt = Deltaopt_cube(1,d,n,s,accuracy)
    % the expected distortion (distortion = (Ls-mean quantization error)^s) is:
        Edist = Edistortion_cube(1, Deltaopt,d,n,s)
    %
    % Here n < 2^30
    % For n vertices sampled among the vertices of [-delta,delta]^d, the optimal delta is:
        Deltaopt = Deltaopt_cube(0,d,n,s,accuracy)
    % the expected distortion (distortion = (Ls-mean quantization error)^s) is:
        Edist = Edistortion_cube(0, Deltaopt,d,n,s)
        
