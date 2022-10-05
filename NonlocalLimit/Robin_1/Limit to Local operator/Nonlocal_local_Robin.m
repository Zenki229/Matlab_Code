clear;clc;close all;

f = @(x) x.*(1-x);
al = 3;
g = [1;1];
 del = [0.01,0.005];
% del = [0.16,0.08,0.04];
%  del = 0.1;
% kernel
s = -0.25; 
kerfun = @kernel;
figure;
hold on 
u = load('sol_local').u;


% scalling
auxfunc = @(x) kerfun(x,s,1).*x.^2;
aux = quadgk(@(x) auxfunc(x),0,1) *2;
Scale = 2/aux;

%FEM
M = 2000;
h = 1/M;
% del = h;
vert = (0:h:1)';
nn = size(vert,1);
t  = [(1:nn-1)' (2:nn)'];
nf = (2:nn-1)';
plot(vert(nf),u(vert(nf)),'LineWidth',1.8,'LineStyle','-','Marker','o','MarkerIndices',1:100:M,'DisplayName','local')
[Mass,SS] = Mass_Stiff_1D(vert,t,nf);
U_gen = u(vert(nf));
for k = 1:length(del)
    delta = del(k);
    g1 = @(x) -g(1)/delta;
    g2 = @(x) g(2)/delta;
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
    S(1:nn1+1,1:nn1+1) = 0;
    S(end-nn2:end,end-nn2:end) = 0;
    for i = 1:nn1+1
        S(i,i) = -sum(S(i,:));
    end
    for i = nn_all-nn2:nn_all
        S(i,i) = -sum(S(i,:));
    end
    Stiff = Scale*S;
    S = Scale*S;
%     M1 = h*diag(ones(nn1+1,1));
%     M2 = h*diag(ones(nn2+1,1));
%     S(1:nn1+1,1:nn1+1) = S(1:nn1+1,1:nn1+1)+1*al*M1/delta;
%     S(end-nn2:end,end-nn2:end)= S(end-nn2:end,end-nn2:end)-1*al*M2/delta;
%     
    M1 = h*diag(ones(nn1,1));
    M2 = h*diag(ones(nn2,1));
    S(1:nn1,1:nn1) = S(1:nn1,1:nn1)+1*al*M1/delta;
    S(end-nn2+1:end,end-nn2+1:end)= S(end-nn2+1:end,end-nn2+1:end)-1*al*M2/delta;
    G1 = load_vector(vert1,t1,(2:nn1+1)',g1);
    G2 = load_vector(vert2,t2,(2:nn2+1)',g2);
    FF = load_vector(vert,t,1:nn,f);
    F=[G1(1:end-1);G1(end)+FF(1);FF(nf);FF(end)+G2(1);G2(2:end)];
    
    if al == 0
        SN = ones(size(S,1)+1,size(S,2));
        SN(1:size(S,1), 1:size(S,2)) = S;
        FN = [F;0];
        U = SN\FN;
    else
        U = S\F;
    end
    if s <=0
        Usol = griddedInterpolant(vert(nf),U(nn1+nf-1));
        plot(0.01:h:1-0.01,Usol(0.01:h:1-0.01),'LineWidth',1.8,'LineStyle','-','DisplayName',['\delta=' num2str(delta)]); 
    else 
    plot(vert(nf),U(nn1+nf-1),'LineWidth',1.8,'LineStyle','-','DisplayName',['\delta=' num2str(delta)]); 
    end
    legend
    err_vec = U(nn1+nf-1)-U_gen;
    err(k) = err_vec'*Mass*err_vec;
    err(k) = sqrt(err(k));
end
vpa(err,4)
Checkratio(err,del(end:-1:1))
sum(Checkratio(err,del(end:-1:1)))/(numel(del)-1)