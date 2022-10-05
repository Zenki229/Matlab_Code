clear;clc;close all;

M = 1000; h = 1/M;
f = load('soce').f;
q = load('pot').q;


vert = (0 :h :1)';
nn = size(vert,1);
t = [(1:nn-1)' (2:nn)'];
nf = (1:nn)';

[Mass,Stiff] = Mass_Stiff_1D(vert,t,nf);

% deal with source operator
F = proj_l2_1D(vert,t,nf,Mass,f);
% plot(vert,F);

%dela with weighted mass mat
Q = Weighted_Mass_1D(vert,t,nf,q);

% Compute terminal value 
U0 = zeros(nn,1);

T = 1; N = 2^10; tau = T/N;
A = Mass + tau*(Stiff+Q);
BCD = zeros(N+1,2);

for i = 1:N
    RHS = tau*(Mass*F)+Mass*U0;
    U1 = conjgrad(A,RHS,1e-10);
    BCD(i+1,:) = [U1(1) U1(end)];
    U0 = U1;
end
plot(vert,U0);
uT = griddedInterpolant(vert,U0);
save('final_value','uT');
Time = 0:tau:T;
bcd{1} = griddedInterpolant(Time,BCD(:,1));
bcd{2} = griddedInterpolant(Time,BCD(:,2));
save('bdry','bcd');