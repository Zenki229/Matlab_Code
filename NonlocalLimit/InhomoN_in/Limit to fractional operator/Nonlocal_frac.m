clear;clc;close all;

f = @(x) -exp(-x);
g = [-1;-exp(-1)];
% f = @(x) 4*pi^2*sin(2*pi .*x);
% g = [2*pi;2*pi]; %Neumann bdry

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
M = 1000;
h = 1/M;
vert = (0:h:1)';
nn = size(vert,1);
t  = [(1:nn-1)' (2:nn)'];
nf = (2:nn-1)';
plot(vert(nf),u(vert(nf)),'LineWidth',1.8,'LineStyle','-','Marker','o','MarkerIndices',1:200:M,'DisplayName','local')
[M,SS] = Mass_Stiff_1D(vert,t,nf);
U_gen = u(vert(nf));
for k = 1:length(del)
    delta = del(k);
    g1 = @(x) -g(1)/delta;
%     g1 = @(x) g(1).*sin(1/delta.*x)./(pi.*x);
    g2 = @(x) g(2)/delta;
%     g2 = @(x) g(2).*sin(1/delta.*(x-1))./(pi.*(x-1));
    row = Stiff_nonlocal_row_free(h, s, delta, kerfun);
    vert1 = (-delta:h:h)';
    nn1 = numel(vert1)-2;
    t1 = [(1:nn1-1)' (2:nn1)'];
    vert2 = (1-h:h:1+delta)';
    nn2 = numel(vert2)-2;
    t2 = [(1:nn2-1)' (2:nn2)'];

    nn_all = nn+nn1+nn2-2;
    r = delta/h;
   
    S = zeros(nn_all,nn_all);
    for i = 1:nn_all
        for j = 1:numel(row)
            if i+j-1>nn_all 
                S(i,i) = S(i,i)+row(j);
            else 
            S(i,i+j-1) = S(i,i+j-1)+row(j);
            end
        end
    end
    for i = 1:nn_all
        for j = 2:numel(row)
            if i-j+1 < 1 
                S(i,i) = S(i,i)+row(j);
            else 
            S(i,i-j+1) = S(i,i-j+1)+row(j);
            end
        end
    end
    G1 = load_vector(vert1,t1,(1:nn1)',g1);
    G2 = load_vector(vert2,t2,(1:nn2)',g2);
%     G1 = g1(vert1(2:end-1));
%     G2 = g2(vert2(2:end-1));
    FF = load_vector(vert,t,1:nn,f);
    F=[G1(1:end-1);G1(end)+FF(1);FF(nf);FF(end)+G2(1);G2(2:end)];

    
    Stiff = ones(nn_all+1,nn_all);
    Stiff(1:nn_all,1:nn_all) = Scale*S;
    ff=[F;0];
    U = Stiff\ff;
    plot(vert(nf),U(nn1+nf-1),'LineWidth',1.8,'LineStyle','-','DisplayName',['\delta=' num2str(delta)]); 
    legend
%     err_vec = U(nn1+nf)-U_gen;
%     err(k) = err_vec'*M*err_vec;
end
% Checkratio(err,del(end:-1:1))