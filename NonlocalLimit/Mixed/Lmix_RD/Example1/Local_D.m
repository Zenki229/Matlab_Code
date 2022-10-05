% Possion Problem with Pure Nuemann Boundary conditions equal to 0
% From compatibility condition int_f = g2-g1
clear;clc;close all;
% f = @(x) 4*pi^2*cos(2*pi.*x); 
% al = 2; be = 1;
% g = [2;2];

% f = @(x) x.*(1-x);
f = @(x) 0.*x;
g = [1;1];


M = 10000;
h = 1/M;

nn = M-1;

F = f(h*(1:nn)');
A = 2*diag(ones(nn,1))- diag(ones(nn-1,1),-1)-diag(ones(nn-1,1),1);

A(1,1) = A(1,1);
A(nn,nn) = A(nn,nn);
A = A/h^2;
F(1) = F(1)+g(1)/h^2;
F(end) = F(end)+g(2)/h^2;

U = A\F;
plot(h*(1:M-1)',U);
u = griddedInterpolant(h*(1:M-1)',U);
save('sol_local','u');
