% Possion Problem with Pure Nuemann Boundary conditions equal to 0
% From compatibility condition int_f = g2-g1
clear;clc;close all;
% f = @(x) -exp(-x);
% g = [-1;-exp(-1)]; %Neumann bdry

f = @(x) x.*(1-x);
% int_f = quadgk(@(x) f(x),0,1);
g = [1/6;0];
% u = @(x) sin(2*pi.*x);

% f = @(x) x.*(1-x);
% g = [0;0];
% int_f = quadgk(@(x) f(x),0,1);

M = 2000;
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
S = zeros(size(Stiff,1)+1,size(Stiff,2));
S(1:end-1,:) = Stiff;
S(end,:) = ones(size(Stiff,2),1);
FF = zeros(size(F,1)+1,1);
FF(1:end-1) = F;
U = S\FF;
plot(vert,U); 
u = griddedInterpolant(vert,U);
save('sol_local','u');
