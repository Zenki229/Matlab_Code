clear;clc;close all;

f = @(x) 4*pi^2*sin(2*pi.*x); 
del = [0.32,0.08,0.02];
% kernel
s = 0.5; 
kerfun = @kernel;
figure;
hold on 
u = load('sol_local').u;


% scalling
auxfunc = @(x) kerfun(x,s,1).*x.^2;
aux = quadgk(@(x) auxfunc(x),0,1) *2;
Scale = 2/aux;

%FEM
M = 100;
h = 1/M;
vert = (0:h:1)';
nn = size(vert,1);
t  = [(1:nn-1)' (2:nn)'];
nf = (2:nn-1)';
plot(vert(nf),u(vert(nf)),'LineWidth',1.8,'LineStyle','-','Marker','o','MarkerIndices',1:4:M,'DisplayName','local')
[Mass,SS] = Mass_Stiff_1D(vert,t,nf);
U_gen = u(vert(nf));
for k = 1:length(del)
    delta = del(k);
    row1 = Stiff_nonlocal_row_free(h, s, delta, kerfun);
    Stiff =toeplitz(row1); 
    r = delta/h;
    FF = load_vector(vert,t,nf,f);
    Stiff = Scale*Stiff;
    U = Stiff\FF;
    plot(vert(nf),U,'LineWidth',1.8,'LineStyle','-','DisplayName',['\delta=' num2str(delta)]); 
    legend
%     err_vec = U(nn1+nf)-U_gen;
%     err(k) = err_vec'*Mass*err_vec;
end
% Checkratio(err,del(end:-1:1))