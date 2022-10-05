clc;clear;close all
%delta=0.2;
u_local = load('u_local').u_local;
Ms=[10,20,40,80,160,320,640,1280];
p=0.6;
kerfun = @kernel;ffun=@(x) x.*(1-x);%.%x.^2.*(1-x).^2;
delta=0.1;

for k=1:length(Ms)
    M=Ms(k);h = 1/M; x=0:h:1; 
   [row]=Stiff_nonlocal(h, p, delta, kerfun);
   Mass=h/6*(sptoeplitz([4,1,zeros(1,M-3)]));
   Stiff{k}=toeplitz(row);   
   f=proj_l2(ffun,x);
    us{k}=Stiff{k}\f;
end
plot(x(2:end-1),us{k})
hold on 
plot(x(2:end-1),u_local(x(2:end-1)))