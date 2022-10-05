clc;clear;close all;

%set data for u = exp(-(x-0.5)^2), consider -10< x< 10

s=0.5;al=2*s;

cns=2^(2*s)*s*gamma(s+1/2)/(pi^(1/2)*gamma(1-s));
ue=@(x) exp(-(x-0.5).^2);%-sqrt(2*pi); 
%since the solution u decays fast, 
%so we may consider finite horizon to compute the source $f$ 
%and the neumann boundary $Nu$.
M=20; h=1/M; % set h, don't need to be too small.

x_int= h:h:1-h;

%% set right-hand-side
f = load('f_frac_50').f;
ff = f(x_int');
plot(x_int, ff);
hold on 
vert = (0:h:1)';
nn = size(vert,1);
t  = [(1:nn-1)' (2:nn)'];
nf = (2:nn-1)';
[Mass,SS] = Mass_Stiff_1D(vert,t,nf);
kerfun = @kernel;
del = [10,20,40,80,160];
for k = 1:length(del)
    delta = del(k);
    row1 = Stiff_nonlocal_row_free(h, s, delta, kerfun);
    Mr = M*delta+1;
    row1 = row1';
    row10= zeros(2*Mr+1,1); row10(1:Mr+1)= row1(end:-1:end-Mr);
    row10(Mr+2:end)=row1(2:end);
    F = zeros(M-1,1);
    for j=1:length(x_int)    
        xx=(x_int(j) - Mr*h):h:(x_int(j) + Mr*h);   
        F(j) =  ue(xx)*row10+h*ue(x_int(j)); 
    end
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
        for j = 1:numel(row1)
            if i+j-1>nn_all 
                S(i,i) = S(i,i)+row1(j);
            else 
            S(i,i+j-1) = S(i,i+j-1)+row1(j);
            end
        end
    end
    for i = 1:nn_all
        for j = 2:numel(row1)
            if i-j+1 < 1 
                S(i,i) = S(i,i)+row1(j);
            else 
            S(i,i-j+1) = S(i,i-j+1)+row1(j);
            end
        end
    end
    S(1:nn1,1:nn1) = 0;
    S(end-nn2+1:end,end-nn2+1:end) = 0;
    for i = 1:nn1
        S(i,i) = -sum(S(i,:));
    end
    for i = nn_all-nn2+1:nn_all
        S(i,i) = -sum(S(i,:));
    end
    S(nn1+1:end-nn2,nn1+1:end-nn2) = S(nn1+1:end-nn2,nn1+1:end-nn2) + h*diag(ones(nn-2,1));
    U = [ue((-delta+h:h:0)');ue(h:h:1-h)';ue(1:h:1+delta-h)'];
    FF = S*U;
    FF = FF(nn1+nf-1);
    FF = FF/h;
    F= F/h;
    errv{k} = (ff)-FF;
    err(k) = norm(errv{k});
    plot(x_int,FF,'DisplayName',num2str(delta));
    legend
end

