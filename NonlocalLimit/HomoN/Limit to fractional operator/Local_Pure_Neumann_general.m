% Possion Problem with Pure Nuemann Boundary conditions equal to 0
% From compatibility condition int_f = g2-g1
clear;clc;close all;
f = @(x) 3+0.*x;
int_f = quadgk(@(x) f(x),0,1);
g = zeros(2,1);
g(2) = int_f-3;
g(1) = 3;



M = 100;
h = 1/M;

vert = (0:h:1)';
vert(end) = 1;
nn = size(vert,1);
t = [(1:nn-1)' (2:nn)'];
nf =  (1:nn)';

[Mass,Stiff] = Mass_Stiff_1D(vert,t,nf);
F = load_vector(vert,t,nf,f);
F(1) = F(1)-g(1);
F(end) = F(end)+g(2);
U = conjgrad(Stiff,F);
plot(vert,U); 

% [Mass,Stiff] = Mass_Stiff_1D(vert,t,nf);
% F = load_vector(vert,t,nf,f);
% c = ones(M+1,1);
% w = Mass*c;
% PT = eye(M+1)- w*c'/dot(w,c);
% 
% U = conjgrad(PT*Stiff*PT',PT*F);
% plot(vert,U); hold on 
% plot(vert,u(vert))