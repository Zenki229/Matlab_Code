clear;clc;close all;

f = @(x) x.*(1-x);
g = [1;1];
del = [0.16,0.08,0.04,0.02,0.01,0.005,0.0025];
% kernel
s = 0.75;
kerfun = @kernel;
figure;
hold on 
u = load('sol_local').u;

g2 = @(x) g(2)+0.*x;
g1 = @(x) g(1)+0.*x;
% scalling
auxfunc = @(x) kerfun(x,s,1).*x.^2;
aux = quadgk(@(x) auxfunc(x),0,1) *2;
Scale = 2/aux;

%FEM
M = 2000;
h = 1/M;
% del = h;
vert = (0:h:1)';
nn = size(vert,1);
t  = [(1:nn-1)' (2:nn)'];
nf = (2:nn-1)';
plot(vert(nf),u(vert(nf)),'LineWidth',1.8,'LineStyle','-','Marker','o','MarkerIndices',1:400:M,'DisplayName','local')
[Mass,SS] = Mass_Stiff_1D(vert,t,nf);
U_gen = u(vert(nf));
for k = 1:length(del)
    delta = del(k);
    row = Stiff_nonlocal_row_free(h, s, delta, kerfun);

    vert1 = (-delta:h:h)';
    nn1 = numel(vert1)-2;
    t1 = [(1:nn1+1)' (2:nn1+2)'];
    vert2 = (1-h:h:1+delta)';
    nn2 = numel(vert2)-2;
    t2 = [(1:nn2+1)' (2:nn2+2)'];

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
    S = Scale *S;
    S(1:nn1,:) = 0;
    S(1:nn1,1:nn1) = diag(ones(nn1,1));
    S(end-nn2+1:end,:) = 0;
    S(end-nn2+1:end,end-nn2+1:end) = diag(ones(nn2,1));
    G1 = g1(vert1(2:nn1+1));
    G2 = g2(vert2(2:nn1+1));
    FF = load_vector(vert,t,1:nn,f);
    F=[G1(1:end-1);G1(end)+FF(1);FF(nf);FF(end)+G2(1);G2(2:end)];
    U = S\F;
    plot(vert(nf),U(nn1+nf-1),'LineWidth',1.8,'LineStyle','-','DisplayName',['\delta=' num2str(delta)]); 
    legend
    err_vec = U(nn1+nf-1)-U_gen;
    err(k) = err_vec'*Mass*err_vec;
    err(k) = sqrt(err(k));
end
vpa(err,4)
Checkratio(err,del(end:-1:1))
sum(Checkratio(err,del(end:-1:1)))/(numel(del)-1)