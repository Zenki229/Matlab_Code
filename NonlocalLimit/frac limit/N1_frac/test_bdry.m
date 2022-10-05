clc;clear;close all;

%set data for u = exp(-(x-0.5)^2), consider -10< x< 10

s=0.5;al=2*s;

cns=2^(2*s)*s*gamma(s+1/2)/(pi^(1/2)*gamma(1-s));
ue=@(x) exp(-(x-0.5).^2);%-sqrt(2*pi); 
%since the solution u decays fast, 
%so we may consider finite horizon to compute the source $f$ 
%and the neumann boundary $Nu$.
M=20; h=1/M; % set h, don't need to be too small.
M0 = 30000;
row=cns*set_stiff_v2(al,M0+2, h)/2;
row0= zeros(2*M0+1,1); row0(1:M0+1)= row(end:-1:end-M0);
row0(M0+2:end)=row(2:end);


%% set right-hand-side
row(1) = row(1)/2;

del = 160;
x_int = h:h:1-h;
    % set boundary data
    delta = del;  xb1 = -delta+h:h:0; Nu1=zeros(length(xb1),1);    
    for i = 1:length(xb1)
           Nu1(i) = (-ue(xb1(i))+ue(x_int)) * row((length(xb1)-i+2):(length(xb1)-i+2)+length(x_int)-1);
           Nu1(i) = Nu1(i)/h;
%            Nu12(i) = -integral(@(x) (-ue(xb1(i))+ ue(x)).*cns.*(x-xb1(i)).^(-1-2*s),0,1,'RelTol',0,'AbsTol',1e-12);
%            Nu12(i) = -quadgk(@(x) (-ue(xb1(i))+ ue(x)).*cns.*(x-xb1(i)).^(-1-2*s),0.0001,1);
    end
%        Nu12(end) = -integral(@(x) (-ue(xb1(end))+ ue(x)).*cns.*(x-xb1(end)).^(-1-2*s),0,1,'RelTol',0,'AbsTol',1e-12);  
%     Nu12(end) = -quadgk(@(x) (-ue(xb1(end))+ ue(x)).*cns.*(x-xb1(end)).^(-1-2*s),0.0000001,1);
%         Nu12(end) =cns*( 1/h*integral(@(y) (ue(0) - ue(y)).*y.^(-2*s),0,h) +  integral(@(y) (ue(0) - ue(y)).* y.^(-1-2*s),h,1));

    g1 = griddedInterpolant(xb1,Nu1);

vert = (0:h:1)';
nn = size(vert,1);
t  = [(1:nn-1)' (2:nn)'];
nf = (2:nn-1)';
[Mass,SS] = Mass_Stiff_1D(vert,t,nf);
kerfun = @kernel;
del = [80,160];
for k = 1:length(del)
    delta = del(k);
    row1 = Stiff_nonlocal_row_free(h, s, delta, kerfun);
    Mr = M*delta+1;
    row1 = row1';
    row10= zeros(2*Mr+1,1); row10(1:Mr+1)= row1(end:-1:end-Mr);
    row10(Mr+2:end)=row1(2:end);
    vert1 = (-delta:h:h)';
    nn1 = numel(vert1)-2; 
    t1 = [(1:nn1+1)' (2:nn1+2)'];
    vert2 = (1-h:h:1+delta)';
    nn2 = numel(vert2)-2;
    t2 = [(1:nn2+1)' (2:nn2+2)'];

    nn_all = nn+nn1+nn2-2;
    r = delta/h;
    
% 
%     xb1 = -delta+h:h:0;
%     xb2 = 1:h:1+delta-h;
%     nn1 = numel(xb1);
%     nn2 = numel(xb2);
%     nn = numel(x_int);
%     nn_all = nn1+nn+nn2;
%     S = zeros(nn1,nn_all);
%     for i = 1:length(xb1)
%         for j = 1:nn
%             if nn1-i+1+j <= numel(row1)
%                 S(i,nn1+j) =  row1(nn1-i+1+j);
%             end
%         end
%         S(i,i) = -sum(S(i,:));
%     end
    vert1 = (-delta-h:h:h)';
    nn1 = numel(vert1)-2;
    t1 = [(1:nn1+1)' (2:nn1+2)'];
    vert2 = (1-h:h:1+delta+h)';
    nn2 = numel(vert2)-2;
    t2 = [(1:nn2+1)' (2:nn2+2)'];
    S = zeros(nn1,nn_all);
    for i = 1:nn_all
        for j = 1:numel(row1)
            if i+j-1>nn_all 
                S(i,i) = S(i,i)+row1(j);
            else 
            S(i,i+j-1) = S(i,i+j-1)+row1(j);
            end
        end
    end
    for i = 1:nn_all
        for j = 2:numel(row1)
            if i-j+1 < 1 
                S(i,i) = S(i,i)+row1(j);
            else 
            S(i,i-j+1) = S(i,i-j+1)+row1(j);
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
    SS = S(1:nn1,:);
    U = [ue((-delta+h:h:0)');ue(h:h:1-h)';ue(1:h:1+delta-h)'];
    G1 = SS*U;
    G1 = G1/h;
    G{k} = griddedInterpolant(-delta+h:h:0,G1);

    err(k) = norm((g1(vert1(2:end-1))-G{k}(vert1(2:end-1))));
    plot(vert1(2:end-1), G1);
    hold on 
    legend
end
plot(vert1(2:end-1),g1(vert1(2:end-1)))
