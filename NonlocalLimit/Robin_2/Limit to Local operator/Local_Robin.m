% Possion Problem with Pure Nuemann Boundary conditions equal to 0
% From compatibility condition int_f = g2-g1
clear;clc;close all;
% f = @(x) 4*pi^2*cos(2*pi.*x); 
% al = 2; be = 1;
% g = [2;2];

f = @(x) x.*(1-x);
al = 3;
g = [1;1];


M = 1000;
h = 1/M;

vert = (0:h:1)';
vert(end) = 1;
nn = size(vert,1);
t = [(1:nn-1)' (2:nn)'];
nf =  (2:nn-1)';

[Mass,S] = Mass_Stiff_1D(vert,t,1:nn);
S = S;
S(1,1) = S(1,1)+al;
S(end,end) = S(end,end)-al;
F = load_vector(vert,t,1:nn,f);

F(1) = F(1)-g(1);
F(end) = F(end)+g(2);
U = S\F;
plot(vert,U);
u = griddedInterpolant(vert,U);
save('sol_local','u');
