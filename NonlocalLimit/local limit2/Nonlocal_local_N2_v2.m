clear;clc;close all;

f = @(x) x.*(1-x);
% int_f = quadgk(@(x) f(x),0,1);
g = [1/6;0];
% f = @(x) cos(2*pi .*x);
% g = [1;-1];
% f = @(x) 4*pi^2*sin(2*pi .*x);
% g = [2*pi;2*pi]; %Neumann bdry

% del = [0.08,0.04,0.02,0.01,0.005];
del = [0.16,0.08,0.04,0.02,0.01];
% del = [0.64,0.32,0.16,0.08,0.04];
del=0.1;
% kernel
s = -1;
kerfun = @kernel;


u = load('sol_local').u;


% scalling
% scalling
auxfunc = @(x) kerfun(x,s,1).*x.^2;
aux = quadgk(@(x) auxfunc(x),0,1) *2;
Scale = 2/aux;
% Scale = 1;
%FEM
M = 20;
h = 1/M;
% del = 5*h;
vert = (0:h:1)';
nn = size(vert,1);
t  = [(1:nn-1)' (2:nn)'];
nf = (2:nn-1)';
figure1 = figure(1);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
set(axes1,'Fontsize',16,'FontWeight','bold');
% xlim(axes1,[-0.16,0.16]);
xlim(axes1,[-0.5,0.5]);
% ylim(axes1,[-0.08,0]);
plot(2.4:0.05:2.5,2.4:0.05:2.5,'HandleVisibility','off');
[Mas,SS] = Mass_Stiff_1D(vert,t,nf);
U_gen = u(vert(nf));
for k = 1:length(del)
    delta = del(k);
    g1 = @(x) -g(1)/delta+0.*x;
%     g1 = @(x) g(1).*sin(1/delta.*x)./(pi.*x);
    g2 = @(x) g(2)/delta+0.*x;
%     g2 = @(x) g(2).*sin(1/delta.*(x-1))./(pi.*(x-1));
    row = Stiff_nonlocal_row_free_v2(M, s, delta, kerfun);
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
    G1 = load_vector(vert1,t1,(2:nn1+1)',g1);
    G2 = load_vector(vert2,t2,(2:nn2+1)',g2);
    FF = load_vector(vert,t,1:nn,f);
    F=[G1(1:end-1);G1(end);FF(nf);G2(1);G2(2:end)];
%     F=[G1;FF;G2];
 

    Stiff = ones(nn_all,nn_all); Stiff(2:end,:) = 0;
    Stiff = Stiff+Scale*S;
%     Stiff = zeros(nn_all,nn_all); Stiff(1,r+1:r+size(M,1)) = 1;
%     Stiff = Stiff+Scale*S;
%     Mass = zeros(nn_all,nn_all);
%     Mass(r+1:r+size(M,1),r+1:r+size(M,1))=M;
%     Stiff = Scale*S+Mass;
    U = Stiff\F;
    kk = round(0.01/h+1e-4);
    interval = linspace(-1,1,2*r+1);
    plot(interval,U((1:2*r+1)'),'LineWidth',1.8,'LineStyle','-','DisplayName',['\delta=' num2str(delta)]); 

    lgd = legend;
    lgd.Location = 'northwest';
    lgd.FontSize = 12;
    lgd.FontWeight = 'bold';
%     err_vec = U(nn1+nf-1)-U_gen;
%     err(k) = err_vec'*M*err_vec;
%     err(k) = sqrt(err(k));
end
% vpa(err,4)
% sum(Checkratio(err,del(end:-1:1)))/(numel(del)-1)