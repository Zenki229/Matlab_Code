clear;clc;close all;

s = 0.25; 

f = @(x) x.*(1-x);
al = 0;
g = [1/6;0];

del = [0.16,0.08,0.04,0.02,0.01,0.005,0.0025];
% del = [0.02,0.01,0.005,0.0025];


kerfun = @kernel;
figure;
hold on 
u = load('sol_local').u;


% scalling
auxfunc = @(x) kerfun(x,s,1).*x.^2;
aux = quadgk(@(x) auxfunc(x),0,1) *2;
Scale = 2/aux;

%Collocation
M = 2000;
h = 1/M;
% m = 2*2.^(3:-1:0);
% del = m*h;
% del = 0.01;
% del = 0.1;
% del = 0.043:-0.001:0.037;
vert = (0:h:1)';
nn = numel(vert);
nf = (2:nn-1)';
plot(vert(nf),u(vert(nf)),'LineWidth',1.8,'LineStyle','-','Marker','o','MarkerIndices',1:100:M,'DisplayName','local')
% [Mass] = Mass_Stiff_1D(vert,t,nf);
U_gen = u(vert(nf));
for k = 1:length(del)
    
    delta = del(k);
%     g1 = @(x) -g(1)/delta;
%     g2 = @(x) g(2)/delta;
    
    row = Stiff_nonlocal_collocation(delta,h, s, kerfun);
    
    r = round(delta/h+0.000001);
    nn_all = r+nn+r;
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
    S(end-r:end,end-r:end) = 0;
    for i = 1:r+1
        S(i,i) = -sum(S(i,:));
    end
    for i = nn_all-r:nn_all
        S(i,i) = -sum(S(i,:));
    end
    Stiff = Scale*S;
    S = Scale*S;
    
    S(1:r+1,1:r+1) = S(1:r+1,1:r+1);
    S(end-r:end,end-r:end)= S(end-r:end,end-r:end);
    G1 = -g(1)/delta*ones(r+1,1);
    FF = f(h*(0:M)');
    G2 = +g(2)/delta*ones(r+1,1);
    F=[G1(1:r);G1(r+1)+FF(1);FF(2:end-1);FF(end)+G2(1);G2(2:r+1)];
    
    S = S(2:end-1,2:end-1);
    F = F(2:end-1);
    SN = ones(size(S,1)+1,size(S,2));
    SN(1:size(S,1), 1:size(S,2)) = S;
    FN = [F;0];
    U = SN\FN;


    if s <=0
        Usol = griddedInterpolant(h*(1:M-1),U(r+1+(1:M-1)'));
        plot(0.01:h:1-0.01,Usol(0.01:h:1-0.01),'LineWidth',1.8,'LineStyle','-','DisplayName',['\delta=' num2str(delta)]); 
    else 
    plot(h*(1:M-1),U(r+(1:M-1)'),'LineWidth',1.8,'LineStyle','-','DisplayName',['\delta=' num2str(delta)]); 
    end
    legend
    err_vec = U(r+(1:M-1)')-U_gen;
    err(k) = norm(err_vec);
end
vpa(err,4)
Checkratio(err,del(end:-1:1))
sum(Checkratio(err,del(end:-1:1)))/(numel(del)-1)