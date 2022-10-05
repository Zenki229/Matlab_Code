clear;clc;close all;
f = @(x) x.*(1-x);
g = [1;1];
M = 2000;
h = 1/M;

vert = (0:h:1)';
vert(end) = 1;
nn = size(vert,1);
t = [(1:nn-1)' (2:nn)'];
nf =  (1:nn)';

[Mass,Stiff] = Mass_Stiff_1D(vert,t,nf);
F = load_vector(vert,t,nf,f);

Stiff(1,:) = 0;
Stiff(1,1) = 1;
Stiff(end,:) = 0;
Stiff(end,end)=1;

F(1) = g(1);
F(end) = g(2);
U = Stiff\F;
plot(vert,U); 
u = griddedInterpolant(vert,U);
save('sol_local','u');
