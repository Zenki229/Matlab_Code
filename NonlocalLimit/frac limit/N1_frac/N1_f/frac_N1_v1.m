clear;clc;close all;

f = load('f_frac_50').f;


% del = [4,8,16,32];
del=[20,40,80,160,320];
% del = del*2;
% del = 1;
s = 0.5;
kerfun = @kernel;
figure1 = figure;
% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
set(axes1,'FontSize',16,'FontWeight','bold');
% ylim(axes1,[0, 2]);
g1 = load('g1_50').g1;
g2 = load('g2_50').g2;

u=@(x) exp(-(x-0.5).^2);

%FEM
M = 20;
h = 1/M;

vert = (0:h:1)';
nn = size(vert,1);
t  = [(1:nn-1)' (2:nn)'];
nf = (2:nn-1)';
plot(vert(nf),u(vert(nf)),'LineWidth',1.8,'LineStyle','-','Marker','o','MarkerIndices',1:2:M-1,'MarkerSize',10,'DisplayName','fractional')
[M,SS] = Mass_Stiff_1D(vert,t,nf);
U_gen = u(vert(nf));
for k = 1:length(del)
    delta = del(k);
    g1a = g1;
    g2a = g2;
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
    S(1:nn1,1:nn1) = 0;
    S(end-nn2+1:end,end-nn2+1:end) = 0;
    for i = 1:nn1
        S(i,i) = -sum(S(i,:));
    end
    for i = nn_all-nn2+1:nn_all
        S(i,i) = -sum(S(i,:));
    end
    S(nn1+1:end-nn2,nn1+1:end-nn2) = S(nn1+1:end-nn2,nn1+1:end-nn2) + h*diag(ones(nn-2,1));
%     G1 = load_vector(vert1,t1,(2:nn1+1)',g1a);
    G1 = h*g1a(vert1(2:nn1+1));
    G2 = h*g2a(vert2(2:nn2+1));
%     G2 = load_vector(vert2,t2,(2:nn2+1)',g2a);
    FF = h*f(vert(nf));
    F=[G1(1:end-1);G1(end);FF;G2(1);G2(2:end)];
%     F = [G1;FF;G2];
    U = S\F;
    plot(vert(nf),U(nn1+nf-1),'LineWidth',1.8,'LineStyle','-','DisplayName',['\delta=' num2str(delta)]);  
    legend('Location','northeast','FontSize',12,'FontWeight','bold')
    err_vec = U(nn1+nf-1)-U_gen;
    Usol{k} = U(nn1+nf-1);
    err(k) = err_vec'*M*err_vec;
    err(k) = sqrt(err(k));
end
vpa(err,4)
Checkratio(err,del)
sum(Checkratio(err,del))/(numel(del)-1)
for k = 1:length(del)-1
    err_vec = Usol{k+1}- Usol{k};
    err2(k) = err_vec'*M*err_vec;
    err2(k) = sqrt(err2(k));
end
vpa(err2,4)
Checkratio(err2,del(1:end-1))
sum(Checkratio(err2,del(1:end-1)))/(numel(del)-2)