clear;clc;close all;
f = @(x) 4*pi^2*sin(2*pi.*x); 
M = 1000;
h = 1/M;

vert = (0:h:1)';
vert(end) = 1;
nn = size(vert,1);
t = [(1:nn-1)' (2:nn)'];
nf =  (2:nn-1)';

[Mass,Stiff] = Mass_Stiff_1D(vert,t,nf);
FF = load_vector(vert,t,nf,f);


U = Stiff\FF;
plot(vert,[0;U;0]); 
u = griddedInterpolant(vert,[0;U;0]);
save('sol_local','u');
