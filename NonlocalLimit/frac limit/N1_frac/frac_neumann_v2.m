clc;clear;close all;

%set data for u = exp(-(x-0.5)^2), consider -10< x< 10

s=0.5;al=2*s;

cns=2^(2*s)*s*gamma(s+1/2)/(pi^(1/2)*gamma(1-s));
ue=@(x) exp(-(x-0.5).^2);%-sqrt(2*pi); 
%since the solution u decays fast, 
%so we may consider finite horizon to compute the source $f$ 
%and the neumann boundary $Nu$.
M=100; h=1/M; % set h, don't need to be too small.

M0=50000; radius = h*M0; %finite horizon to set f and Nu
                   
row=cns*set_stiff_v2(al,M0+2, h)/2;
row0= zeros(2*M0+1,1); row0(1:M0+1)= row(end:-1:end-M0);
row0(M0+2:end)=row(2:end);

x_int= h:h:1-h;

%% set right-hand-side
% ff = zeros(M-1,1);
% for j=1:length(x_int)    
%     xx=(x_int(j) - M0*h):h:(x_int(j) + M0*h);   
%     ff(j) =  ue(xx)*row0+h*ue(x_int(j));
%     
% end
% 
% ff = ff/h;
del=320;
% del = del*2;



%% set boundary data
row(1) = row(1)/2;


    % set boundary data
    delta = del;  xb1 = -delta+h:h:0; Nu1=zeros(length(xb1),1);    
    for i = 1:length(xb1)
           Nu1(i) = (-ue(xb1(i))+ue(x_int)) * row((length(xb1)-i+2):(length(xb1)-i+2)+length(x_int)-1);
           Nu1(i) = Nu1(i);
           Nu12(i) = -integral(@(x) (-ue(xb1(i))+ ue(x)).*cns.*(x-xb1(i)).^(-1-2*s),0,1,'RelTol',0,'AbsTol',1e-12);
%            Nu12(i) = -quadgk(@(x) (-ue(xb1(i))+ ue(x)).*cns.*(x-xb1(i)).^(-1-2*s),0.0001,1);
    end
%        Nu12(end) = -integral(@(x) (-ue(xb1(end))+ ue(x)).*cns.*(x-xb1(end)).^(-1-2*s),0,1,'RelTol',0,'AbsTol',1e-12);  
%     Nu12(end) = -quadgk(@(x) (-ue(xb1(end))+ ue(x)).*cns.*(x-xb1(end)).^(-1-2*s),0.0000001,1);
        Nu12(end) =cns*( 1/h*integral(@(y) (ue(0) - ue(y)).*y.^(-2*s),0,h) +  integral(@(y) (ue(0) - ue(y)).* y.^(-1-2*s),h,1));

    g1 = griddedInterpolant(xb1,Nu12);
    Nu22 = Nu12(end:-1:1);
    xb2 = 1:h:1+delta-h;
    g2 = griddedInterpolant(xb2,Nu22);
% for k=1:length(del)    
%     % set boundary data
%     delta = del(k);  xb2 = 1:h:1+delta-h; Nu2=zeros(length(xb2),1);  
%     for i = 1:length(xb2)
%         j = length(xb2)-i+1;
%             Nu2(i) = (-ue(xb2(i))+ue(x_int(end:-1:1))) * row((length(xb2)-j+2):(length(xb2)-j+2)+length(ff)-1);
%             Nu2(i) = Nu2(i);
%     end
%     g2{k} = griddedInterpolant(xb2,Nu2);
% end
f = griddedInterpolant(x_int,ff);
save('f_frac_50','f');
% save('g1_16','g1');
% save('g2_16','g2');
save('g1_50','g1');
save('g2_50','g2');



%h=1:100


 



% x=-10:0.1:10;
% 
% plot(x,ue(x))