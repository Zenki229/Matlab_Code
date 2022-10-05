function [F] = proj_H1_1D(vert,t,nf,Stiff,f)
[node,W] = gauleg(0,1,20);
nn = size(vert,1);
nt = size(t,1);
h = (1/nn);
df = @(x) (f(x+h)-f(x-h))./(2*h);
b = zeros(nn,1);
for i = 1:nt
    x0 = vert(t(i,1));
    x1 = vert(t(i,2));
    fy(1) =  -1/(x1-x0);
    fy(2) =  1/(x1-x0);
    aux = zeros(2,1);
    
    for j = 1:2
        u = @(x) df((x1-x0).*x+x0);
        aux(j) =  (x1-x0)*fy(j)*dot(u(node),W);
%         aux(j) = quad(@(x) fy(j).*df(x),x0,x1);
    end
    b(t(i,:)) = b(t(i,:))+aux;
end
b = b(nf);
 F = conjgrad(Stiff,b);
