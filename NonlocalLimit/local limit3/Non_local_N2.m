clear;clc;close all;

f = @(x) x.*(1-x);
g = [1/6;0];
% del = [0.08,0.04,0.02,0.01,0.005];
% del = [0.16,0.08,0.04];
del =0.1;
s = -1;
kerfun = @kernel;
u = load('sol_local').u;
auxfunc = @(x) kerfun(x,s,1).*x.^2;
aux = quadgk(@(x) auxfunc(x),0,1) *2;
Scale = 2/aux;
M = 20;
h = 1/M;
vert = (0:h:1)';
nn = size(vert,1);
t  = [(1:nn-1)' (2:nn)'];
nf = (2:nn-1)';
figure1 = figure(1);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
set(axes1,'Fontsize',16,'FontWeight','bold');
plot(2.4:0.05:2.5,2.4:0.05:2.5,'HandleVisibility','off');
% plot(vert(nf),u(vert(nf)),'LineWidth',1.8,'LineStyle','-','Marker','o','MarkerIndices',1:100:M,'DisplayName','local')
xlim(axes1,[-0.5,1]);
% ylim(axes1,[-0.06,0.04]);
[Mas,SS] = Mass_Stiff_1D(vert,t,nf);
U_gen = u(vert(nf));
normu = U_gen'*Mas*U_gen;normu=sqrt(normu);
global x_int xb1 xb2;
for k = 1:length(del)
    delta = del(k);
    r=floor(delta/h+1e-5);
    g1 = @(x) -g(1)/delta+0.*x;
    g2 = @(x) g(2)/delta+0.*x;
    row = Stiff_nonlocal_row_free_v2(M, s, delta, kerfun);
    xb1 = (-delta:h:0)';
    xb2 = (1:h:1+delta)';
    x_int = (0:h:1)';
    nn1 = numel(xb1);
    nn2 = numel(xb2);
    nn = numel(x_int);
    nn_all = nn1+nn2+nn;
    G1 = setRHS(g1,xb1);
    G11 = h*g1(xb1);
    G2 = h*g2(xb2);
    FF = load_vector(vert,t,(1:size(vert,1))',f);
    F = [G1;FF;G2];
%     F = [G1(1:end-1);FF(1)+G1(end);FF(nf);FF(end)+G2(1);G2(2:end)];
    S = Stiff_nonlocal_mat_v2(row,h,delta);
    Stiff = ones(nn_all,nn_all); Stiff(2:end,:) = 0;
    Stiff = Stiff+Scale*S;
    U = Stiff\F;
    interval = linspace(-1,1,2*r+1);
    plot(interval,U((1:2*r+1)'),'LineWidth',1.8,'LineStyle','-','DisplayName',['\delta=' num2str(delta)]); 
%     plot(x_int,U((nn1+1:nn1+nn)'),'LineWidth',1.8,'LineStyle','-','DisplayName',['\delta=' num2str(delta)]); 
    lgd = legend;
    lgd.Location = 'northwest';
    lgd.FontSize = 12;
    lgd.FontWeight = 'bold';
%     err_vec = U(nn1+nf-1)-U_gen;
%     err(k) = err_vec'*Mas*err_vec;
%     err(k) = sqrt(err(k));
end
% vpa(err,4)
% sum(Checkratio(err,del(end:-1:1)))/(numel(del)-1)