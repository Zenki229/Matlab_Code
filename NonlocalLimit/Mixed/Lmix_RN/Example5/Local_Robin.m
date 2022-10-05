% Possion Problem with Pure Nuemann Boundary conditions equal to 0
% From compatibility condition int_f = g2-g1
clear;clc;close all;
% f = @(x) 4*pi^2*cos(2*pi.*x); 
% al = 2; be = 1;
% g = [2;2];

f = @(x) x.*(1-x);
al = 2;
g = [1;1/2];


M = 4000;
h = 1/M;

vert = (0:h:1)';
nn = numel(vert);

F = f(vert);
A = 2*diag(ones(nn,1))- diag(ones(nn-1,1),-1)-diag(ones(nn-1,1),1);

A(1,1) = A(1,1)-(1-al*h);
A(nn,nn) = A(nn,nn)-(1+al*h);
A = A/h^2;
F(1) = F(1)- g(1)/h;
F(end) = F(end) + g(2)/h;

U = A\F;
plot(vert,U);
u = griddedInterpolant(vert,U);
save('sol_local','u');
