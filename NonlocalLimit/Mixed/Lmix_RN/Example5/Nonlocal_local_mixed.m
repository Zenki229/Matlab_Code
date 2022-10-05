clear;clc;close all;

s = -1; 
%-Lu + au = f
% N(u) = -3/delta, (-delta,-delta/2)
% u = 0 ,(-delta/2,0)
% N(u) = 3/delta (1,1+delta)


% f = @(x) x.*(1-x);
f = @(x) 0.*x;
a = 1;
g = [-1/3;0;1/3];



% del = [0.02,0.01];
del = [0.16,0.08,0.04,0.02,0.01];
% del = [0.16,0.04,0.01,0.0025];
% del = del/2;
kerfun = @kernel;


figure1 = figure;
% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
set(axes1,'FontSize',16,'FontWeight','bold');
% ylim(axes1,[0, 2]);

u = load('sol_local').u;
uN = load('sol_localN').u;

% scalling
auxfunc = @(x) kerfun(x,s,1).*x.^2;
aux = quadgk(@(x) auxfunc(x),0,1) *2;
Scale = 2/aux;

%Collocation
M = 10000;
h = 1/M;
% del = 0.01;
plot(h*(1:M-1),u(h*(1:M-1)),'LineWidth',1.8,'LineStyle','-','Marker','o','MarkerIndices',1:1000:M,'DisplayName','local')
% plot(h*(1:M-1),uN(h*(1:M-1)),'LineWidth',1.8,'LineStyle','-','Marker','+','MarkerIndices',1:1000:M,'DisplayName','localN')

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
    S = Scale*S;
    Stiff = Scale*S;

%     r1 = floor((r+1)/2);
%     r2 = r+1-r1;
    r2 = 3;
    r1 = r+1-r2;
    S(1:r+1,1:r+1) = 0;
     for i = 1:r1
        S(i,i) = -sum(S(i,:));
    end
    for i = r1+1:r+1
        S(i,:) = 0;
        S(i,i) = 1;
    end
    S(nn_all-r:nn_all,nn_all-r:nn_all) = 0;
    for i = nn_all-r:nn_all
       S(i,i) = -sum(S(i,:));
    end
    A = a * diag(ones(nn_all,1));
    G1 = g(1)/delta * ones(r1,1); 
    G2 = g(2)*ones(r2,1);
    FF = f(h*(0:M)');
    G3 = g(3)/delta * ones(r+1,1);
    F=[G1;G2(1:end-1);G2(end)+FF(1);FF(2:end-1);FF(end)+G3(1);G3(2:r+1)];
    SA = S+A;
    U = SA(2:end-1,2:end-1)\F(2:end-1);
    plot(h*(1:M-1),U(r+(1:M-1)'),'LineWidth',1.8,'LineStyle','-','DisplayName',['\delta=' num2str(delta)]); 
    legend('location','northwest')
    err_vec = U(r+(1:M-1)')-U_gen;
    Usol{k} = U(r+(1:M-1)');
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