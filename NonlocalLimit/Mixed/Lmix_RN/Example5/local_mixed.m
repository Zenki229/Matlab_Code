clear;clc;close all;

% f = @(x) x.*(1-x);
f = @(x) 0.*x;
al = 0;
g = [0;1/3];


M = 10000;
h = 1/M;


nn = M-1;

F = f(h*(1:nn)');
A = 2*diag(ones(nn+1,1))- diag(ones(nn,1),-1)-diag(ones(nn,1),1);
A(end,end) = A(end,end)-1;
A = A/h^2;
A(1,:) = 0;
A(1,1) = 1;
Da = diag(ones(nn+1,1));
Da(1,:) = 0;


F(end) = F(end) + g(2)/h;
FF =[g(1);F];
U = (A+Da)\FF;
plot(h*(1:nn)',U(2:end));
u = griddedInterpolant(h*(1:nn)',U(2:end));
save('sol_local','u');
