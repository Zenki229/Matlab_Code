% Possion Problem with Pure Nuemann Boundary conditions equal to 0
% From compatibility condition int_f = g2-g1
clear;clc;close all;
f = @(x) x.*(1-x);


M = 1000;
h = 1/M;

vert = (0:h:1)';
vert(end) = 1;
nn = size(vert,1);
t = [(1:nn-1)' (2:nn)'];
nf =  (2:nn-1)';

[Mass,Stiff] = Mass_Stiff_1D(vert,t,nf);
F = load_vector(vert,t,nf,f);
U = conjgrad(Stiff,F);
plot(vert(nf),U)
u_local = griddedInterpolant(vert(nf),U);
save('u_local','u_local');

% [Mass,Stiff] = Mass_Stiff_1D(vert,t,nf);
% F = load_vector(vert,t,nf,f);
% c = ones(M+1,1);
% w = Mass*c;
% PT = eye(M+1)- w*c'/dot(w,c);
% 
% U = conjgrad(PT*Stiff*PT',PT*F);
% plot(vert,U); hold on 
% plot(vert,u(vert))