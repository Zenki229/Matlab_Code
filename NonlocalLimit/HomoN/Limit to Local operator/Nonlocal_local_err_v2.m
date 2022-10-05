clear;clc;close all;

f = @(x) cos(2*pi.*x);
del = [0.64,0.32,0.16,0.08,0.04,0.02];
% kernel
s = 0.3;
kerfun = @kernel;
u = @(x) 1/(4*pi^2).*cos(2*pi.*x);


% scalling
auxfunc = @(x) kerfun(x,s,1).*x.^2;
aux = quadgk(@(x) auxfunc(x),0,1) *2;
Scale = 2/aux;

%FEM
M = 2000;
h = 1/M;
vert = (0:h:1)';
nn = size(vert,1);
t  = [(1:nn-1)' (2:nn)'];
nf = (2:nn-1)';
[M,SS] = Mass_Stiff_1D(vert,t,nf);
Ug = u(vert(nf));
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
    S(1:nn1+1,1:nn1+1) = 0;
    S(end-nn2:end,end-nn2:end) = 0;
    for i = 1:nn1+1
        S(i,i) = -sum(S(i,:));
    end
    for i = nn_all-nn2:nn_all
        S(i,i) = -sum(S(i,:));
    end

    FF = load_vector(vert,t,nf,f);
    F = zeros(nn_all,1);
    F(nn1+nf) = FF;
    Stiff = zeros(nn_all+1,nn_all);
    Stiff(1:nn_all,1:nn_all) = S;
    Stiff(end,:) = ones(1,nn_all);
    ff = zeros(nn_all+1,1);
    ff(1:nn_all) = F;
    U = Stiff\ff;
    U = U/Scale;
    err_vec = U(nn1+nf)-Ug;
    err(k) = err_vec'*M*err_vec;
    err(k) = sqrt(err(k));
end
Checkratio(err,del(end:-1:1))