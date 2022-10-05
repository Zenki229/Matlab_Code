clear;clc;close all;

s = -1; 

f = @(x) 0.*x;
al = 2;
g = [-1/3;1;1];




del = [0.16,0.08,0.04,0.02,0.01];
% del = [0.02,0.01,0.005,0.0025];
% del = [0.16,0.08,0.04];

% del = del/2;

kerfun = @kernel;




u = load('sol_local').u;


% scalling
auxfunc = @(x) kerfun(x,s,1).*x.^2;
aux = quadgk(@(x) auxfunc(x),0,1) *2;
Scale = 2/aux;

%Collocation
M = 10000;
h = 1/M;
figure1 = figure(1);
axes1 = axes('Parent',figure1);
hold(axes1,'on');
set(axes1,'Fontsize',16,'FontWeight','bold');
xlim(axes1,[-1,2]);
ylim(axes1,[0.998,1.03]);
% interval2 = linspace(0,2,M-1);
% plot(interval2,u(h*(1:M-1)),'LineWidth',1.8,'LineStyle','-','Marker','o','MarkerIndices',round(linspace(1,M-1,10)),'DisplayName','local')
plot(2.4:0.05:2.5,2.4:0.05:2.5,'HandleVisibility','off');
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

    r1 = ceil((r+1)/2);
    r2 = r+1-r1;
    S(1:r+1,1:r+1) = 0;
    for i = 1:r1
        S(i,i) = -sum(S(i,:));
    end
    for i = r1+1:r+1
        S(i,:) = 0;
        S(i,i) = 1;
    end
    for i = nn_all-r:nn_all
        S(i,:) = 0;
        S(i,i) = 1;
    end
    
    
    G1 = -g(1)/delta*ones(r1,1);
    G2 = g(2) * ones(r2,1); 
    FF = f(h*(0:M)');
    G3 = +g(3) * ones(r+1,1);
    F=[G1;G2(1:end-1);G2(end)+FF(1);FF(2:end-1);FF(end)+G3(1);G3(2:r+1)];
    U = S(2:end-1,2:end-1)\F(2:end-1);
    interval = linspace(-1,2,3*r+1);
    plot(interval,U((1:r+1+2*r)'),'LineWidth',1.8,'LineStyle','-','DisplayName',['\delta=' num2str(delta)]); 

    lgd = legend;
    lgd.Location = 'northeast';
    lgd.FontSize = 12;
    lgd.FontWeight = 'bold';
%     err_vec = U(r+1+(1:M-1)')-U_gen;
%     Usol{k} = U(r+1+(1:M-1)');
%     err(k) = norm(err_vec);
%     err(k) = err(k)*sqrt(h);
end
% vpa(err,4)
% Checkratio(err,del(end:-1:1))
% sum(Checkratio(err,del(end:-1:1)))/(numel(del)-1)
% for k = 1:length(del)-1
%     err_vec = Usol{k+1}- Usol{k};
%     err2(k) = norm(err_vec)*sqrt(h);
% end
% vpa(err2,4)
% Checkratio(err2,del(end:-1:2))
% sum(Checkratio(err2,del(end:-1:2)))/(numel(del)-2)