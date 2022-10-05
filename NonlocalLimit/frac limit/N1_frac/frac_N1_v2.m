clear;clc;close all;

f = load('f_frac').f;


% del = [4,8,16,32];
del=[10,20,40,80,160];
% del = del*2;
% del = 2;
s = 0.25;
kerfun = @kernel;
figure1 = figure;
% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
set(axes1,'FontSize',16,'FontWeight','bold');
% ylim(axes1,[0, 2]);
g1 = load('g1').g1;
g2 = load('g2').g2;

u=@(x) exp(-(x-0.5).^2);

%FEM
M = 50;
h = 1/M;

vert = (0:h:1)';
nn = size(vert,1);
t  = [(1:nn-1)' (2:nn)'];
nf = (2:nn-1)';
plot(vert(nf),u(vert(nf)),'LineWidth',1.8,'LineStyle','-','Marker','o','MarkerIndices',[1:10:M-1,M-1],'MarkerSize',10,'DisplayName','fractional')
[M,SS] = Mass_Stiff_1D(vert,t,nf);
U_gen = u(vert(nf));
for k = 1:length(del)
    delta = del(k);
    g1a = g1{k};
    g2a = g2{k};
    row = Stiff_nonlocal_row_free(h, s, delta, kerfun);
    xb1 = -delta+h:h:0;
    xb2 = 1:h:1+delta-h;
    x_int = h:h:1-h;
    nn1 = numel(xb1);
    nn2 = numel(xb2);
    nn = numel(x_int);
    S = Stiff_nonlocal_mat(row,h,delta);
    S(nn1+(1:nn),nn1+(1:nn)) = S(nn1+(1:nn),nn1+(1:nn))+h*diag(ones(nn,1));
    
    G1 = h*g1a(xb1');
    G2 = h*g2a(xb2');
    FF = h*f(x_int');
    F = [G1;FF;G2];
    
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