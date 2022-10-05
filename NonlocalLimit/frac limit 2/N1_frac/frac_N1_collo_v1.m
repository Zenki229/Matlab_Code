clear;clc;close all;

f = load('f_frac').f;


del = [1,2,4,8,16];
% del = 1;
s = 0.75;
kerfun = @kernel;
figure1 = figure;
% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
set(axes1,'FontSize',16,'FontWeight','bold');
% ylim(axes1,[0, 2]);
g1 = load('g1').g1;
g2 = load('g2').g2;

u=@(x) exp(-(x-0.5).^2/2);
%Collocation
M = 200;
h = 1/M;
% del = 0.1;
plot(h*(1:M-1),u(h*(1:M-1)),'LineWidth',1.8,'LineStyle','-','Marker','o','MarkerIndices',1:1000:M,'DisplayName','local')

U_gen = u(h*(1:M-1)');
for k = 1:length(del)
    
    delta = del(k);
    row = Stiff_nonlocal_collocation(delta,h, s, kerfun);
    
    r = round(delta/h+0.000001);
    nn_all = r+M+1+r;
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




    S(1:r+1,1:r+1) = 0;
    for i = 1:r+1
        S(i,i) = -sum(S(i,:));
    end
    S(nn_all-r:nn_all,nn_all-r:nn_all) = 0;
    for i = nn_all-r:nn_all
       S(i,i) = -sum(S(i,:));
    end
    Stiff = S;
    G1 = g1{k}(h*(-r:1:0)');
    FF = f(h*(0:M)');
    G2 = g2{k}(1+(0:1:r)'*h);
    F=[G1(1:end-1);G1(end)+FF(1);FF(2:end-1);FF(end)+G2(1);G2(2:r+1)];
    S(r+2:end-r-1,r+2:end-r-1) = S(r+2:end-r-1,r+2:end-r-1) + diag(ones(M-1,1));
%     S = S + diag(ones(nn_all,1));
    U = S(1:end-1,1:end-1)\F(1:end-1);
    plot(h*(1:M-1),U(r+(1:M-1)'),'LineWidth',1.8,'LineStyle','-','DisplayName',['\delta=' num2str(delta)]); 
    legend('location','northwest')
    err_vec = U(r+(1:M-1)')-U_gen;
    Usol{k} = U(r+1+(1:M-1)');
    err(k) = norm(err_vec);
    err(k) = err(k)*sqrt(h);
end
vpa(err,4)
Checkratio(err,del(end:-1:1))
sum(Checkratio(err,del(end:-1:1)))/(numel(del)-1)
for k = 1:length(del)-1
    err_vec = Usol{k+1}- Usol{k};
    err2(k) = norm(err_vec)*sqrt(h);
end
vpa(err2,4)
Checkratio(err2,del(end:-1:2))
sum(Checkratio(err2,del(end:-1:2)))/(numel(del)-2)