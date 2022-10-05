clear;clc;close all;

f = @(x) x.*(1-x);
del = [100];
% kernel
s = 0.7;
kerfun = @kernel;

figure;
hold on


%FEM
M = 100;
h = 1/M;
vert = (0:h:1)';
nn = size(vert,1);
t  = [(1:nn-1)' (2:nn)'];
nf = (2:nn-1)';
for k = 1:length(del)
    delta = del(k);
    row = Stiff_nonlocal_row_free(h, s, delta, kerfun);
    nn1 = numel(-delta:h:0);
    nn2 = numel(1:h:1+delta);
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
    FF = load_vector(vert,t,nf,f);
    F = zeros(nn_all,1);
    F(nn1+nf) = FF;
    Stiff = zeros(nn_all+1,nn_all);
    Stiff(1:nn_all,1:nn_all) = S;
    Stiff(end,:) = ones(1,nn_all);
    ff = zeros(nn_all+1,1);
    ff(1:nn_all) = F;
    U = Stiff\ff;
    plot(vert(nf), U(nn1+nf)) 
end
u_sol = griddedInterpolant(vert(nf), U(nn1+nf));