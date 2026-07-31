function [ Energy, Physical_energy ] = Energy_Xn( Xn, theta,p,kernel,TENSORISED)
% function [ Energy, Physical_energy ] = Energy_Xn( Xn, theta,p,kernel,TENSORISED)
%__________________________________________________________________________
% For X (d*n) an n-point design in [0,1]^d, returns Energy and Physical_energy,
% two 1*n vectors corresponding to the energy for the empirical measure associated with Xk
% and the physical energy for Xk, for k=1 (gives NaN),2,...,n
% The kernel used is specified by theta,p,kernel,TENSORISED (see calcR_switch.m)
% Uses calcR_switch.m
%-----
% Author: L. Pronzato, 2022 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE.

[d,n]=size(Xn);
Energy=NaN(1,n);
Physical_energy=NaN(1,n);

Kn=calcR_switch(Xn,theta,p,kernel,TENSORISED);

energy=Kn(1,1); Energy(1)=energy;
physical_energy=0;
for k=2:n
    energy_update=sum(Kn(1:k,k)); energy=energy+energy_update;
    Energy(k)=2*energy/k^2;
    physical_energy_update=sum(Kn(1:k-1,k)); 
    physical_energy=physical_energy+physical_energy_update;
    Physical_energy(k)=2*physical_energy/(k*(k-1));
end
end
