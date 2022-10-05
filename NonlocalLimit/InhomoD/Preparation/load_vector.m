function [F] = load_vector(vert,t,nf,f)
[node,W] = gauleg(0,1,20);
nn = size(vert,1);
nt = size(t,1);

b = zeros(nn,1);
for i = 1:nt
    x0 = vert(t(i,1));
    x1 = vert(t(i,2));
    if x0 == x1
        continue;
    end
    fy{1} = @(x) (x1-x)./(x1-x0);
    fy{2} = @(x) (x-x0)./(x1-x0);
    aux = zeros(2,1);
    
    for j = 1:2
        u = @(x) f((x1-x0).*x+x0).*fy{j}((x1-x0).*x+x0);
        aux(j) =  (x1-x0)*dot(u(node),W);
%         aux(j) = quad(@(x) f(x).*fy{j}(x),x0,x1);
    end
    b(t(i,:)) = b(t(i,:))+aux;
end
b = b(nf);
F = b;
