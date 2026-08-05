function [Gen_letters] = generators(Gen)
% function [Gen_letters] = generators(Gen)
%__________________________________________________________________________
% Gen is a a (m*(d-m)) matrix of generators created by greedy_generators.m or
% SA_exchange_generators.m, for a 2^(d-m) fractional factorial design up to 
% d-m=52. In each row from 1 to m, the first elements of Gen(i,:) equal to
% zero are removed and the next ones indicate the generator for that variable.
% Returns the list of generators using the notation with Latin alphabet: 
%   a for 1, b for 2, z for 26,..., A for 27,..., Z for 52. 
%-----
% Author: L. Pronzato, 2020 <pronzato@i3s.unice.fr>
% Permission is granted to use, modify and redistribute this software for 
% non-commercial research and educational purposes only. Commercial use 
% requires prior written permission from the author.
% This MATLAB function is provided in the hope that it will prove useful, 
% but WITHOUT ANY WARRANTY OF FITNESS FOR A PARTICULAR PURPOSE.

Alphabet='abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUV';
[m,~]=size(Gen);
Gen_letters=cell(m,1);
for i=1:m
    gen_i=Gen(i,:);
    gen_i=gen_i(gen_i>0);
    Gen_letters{i}=Alphabet(gen_i);
end 
end