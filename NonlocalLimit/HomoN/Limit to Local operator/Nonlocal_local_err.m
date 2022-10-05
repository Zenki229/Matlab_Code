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
M = 1000;
h = 1/M;
vert = (0:h:1)';
nn = size(vert,1);
t  = [(1:nn-1)' (2:nn)'];
nf = (2:nn-1)';
[M,SS] = Mass_Stiff_1D(vert,t,nf);
Ug = u(vert(nf));
for k = 1:length(del)
    delta = del(k);
    row1 = Stiff_nonlocal_row_free(h, s, delta, kerfun);
    row = zeros(nn-2,1);
    row(1:numel(row1)) = row1;
    Stiff =toeplitz(row);   
    nn1 = numel(-delta:h:0);
    nn2 = numel(1:h:1+delta);
    nn_all = nn+nn1+nn2;
    r = delta/h;
    S = zeros(nn_all,nn_all);
    S(nn1+1:nn1+nn-2,nn1+1:nn1+nn-2) = Stiff;
    for i = 1:nn1
        S(i,i:i+numel(nn1+i:end)-1) = S(nn1+i,nn1+i:end);
    end

    for i = 1:nn1+r+1
        for j = 1:r+1
            if i - j < 1
                S(i,i) = S(i,i)+S(i,i+j);
                continue;
            end
            S(i,i-j) = S(i,i+j);
        end
    end
    for i = nn1+nn-1:nn_all
        S(i,i-r-1:i) = S(i-1,i-1-r-1:i-1);
    end
    for i = nn1+1:nn_all
        for j = 1:r+1
            if i+j>nn_all
                S(i,i) = S(i,i)+S(i,i-j);
                continue;
            end
            S(i,i+j)=S(i,i-j);
        end
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