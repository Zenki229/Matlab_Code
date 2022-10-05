clear;clc;close all;

% f = @(x) -exp(-x);
% g = [-1;-exp(-1)];

f = @(x) cos(2*pi.*x);
u = @(x) 1/(4*pi^2).*cos(2*pi.*x);
% f = @(x) 2.*x-1;
% u = @(x)1/2.*x.^2-1/3.*x.^3-1/12;
% int_f = quadgk(@(x) f(x),0,1);
g = [0;0];

del = [0.08,0.04,0.02,0.01,0.005];
% del = [0.16,0.08,0.04];
% kernel
s = 0.25;
kerfun = @kernel;
figure1 = figure;
% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');
set(axes1,'FontSize',16,'FontWeight','bold');
% ylim(axes1,[0, 2]);



% scalling
auxfunc = @(x) kerfun(x,s,1).*x.^2;
aux = quadgk(@(x) auxfunc(x),0,1) *2;
Scale = 2/aux;
% Scale = 1;
%FEM
M = 2000;
h = 1/M;
% del = 0.1;
vert = linspace(0,1,M+1); vert = vert';
nn = size(vert,1);
t  = [(1:nn-1)' (2:nn)'];
nf = (1:nn)';
plot(vert,u(vert),'LineWidth',1.8,'LineStyle','-','Marker','o','MarkerIndices',1:200:M,'DisplayName','local')
[Mass,SS] = Mass_Stiff_1D(vert,t,nf);
U_gen = u(vert(nf));
for k = 1:length(del)
    delta = del(k);
    row = Stiff_nonlocal_row_free(h, s, delta, kerfun);
    vert1 = linspace(-delta,0,round(M*delta)+1);vert1 = vert1';
    nn1 = numel(vert1);
    t1 = [(1:nn1-1)' (2:nn1)'];
    vert2 = vert1+1+delta;
    nn2 = numel(vert2);
    t2 = [(1:nn2-1)' (2:nn2)'];

    nn_all = nn+nn1+nn2;
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
    
    FF = load_vector(vert,t,1:nn,f);
    F=[zeros(nn1,1);FF;zeros(nn2,1)];

    
    Stiff = ones(nn_all+1,nn_all);
    Stiff(1:nn_all,1:nn_all) = Scale*S;
    ff=[F;0];
    U = Stiff\ff;
    if s <=0
        Usol = griddedInterpolant(vert(nf),U(nn1+nf-1));
        plot(0.01:h:1-0.01,Usol(0.01:h:1-0.01),'LineWidth',1.8,'LineStyle','-','DisplayName',['\delta=' num2str(delta)]); 
    else 
        plot(vert(nf),U(nn1+nf),'LineWidth',1.8,'LineStyle','-','DisplayName',['\delta=' num2str(delta)]); 
    end    
    legend('Location','southeast')
    err_vec = U(nn1+nf)-U_gen;
    err(k) = err_vec'*M*err_vec;
    err(k) = sqrt(err(k));
end
vpa(err,4)
Checkratio(err,del(end:-1:1))
sum(Checkratio(err,del(end:-1:1)))/(numel(del)-1)