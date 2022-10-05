clc;clear;close all;
%delta=0.2;
u_local = load('u_local').u_local;
 
Ms=10;
p=0.2;
kerfun = @kernel;ffun=@(x) x.*(1-x);%.%x.^2.*(1-x).^2;
% scalling
[node,W] = gauleg(0,1,20);
auxfunc = @(x) kerfun(x,p,1).*x.^2;
% aux = 2*dot(auxfunc(node),W)*2;
aux = quadgk(@(x) auxfunc(x),0,1) *2;
Scale = 2/aux;

kerf = @(x,s,delta) Scale * kernel(x,s,delta);
delta=[0.2];
delta = delta(end:-1:1);
M=Ms;h = 1/M; x=0:h:1; 
plot(x(2:end-1),u_local(x(2:end-1)))
hold on 
for k=1:length(delta)   
   [row]=Stiff_nonlocal(h, p, delta(k), kerfun);
   Mass=h/6*(sptoeplitz([4,1,zeros(1,M-3)]));
   Stiff{k}=toeplitz(row);
   Stiff{k}=Stiff{k}*Scale;
   f=proj_l2(ffun,x);
   us{k}=Stiff{k}\f;
   plot(x(2:end-1),us{k})
end
