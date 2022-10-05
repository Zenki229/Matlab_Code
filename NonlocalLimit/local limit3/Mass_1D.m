function [Mass]= Mass_1D(vert)

nn = size(vert,1);
t = [(1:nn-1)' (2:nn)'];
nt = size(t,1);
nf = (1:nn)';
Mass = zeros(nn,nn);
[node,W] = gauleg(0,1,20);
for i = 1:nt
    x0 = vert(t(i,1));
    x1 = vert(t(i,2));
    if x0 == x1
        continue;
    end
    fy{1} = @(x) (x-x0)./(x1-x0);
    fy{2} = @(x) (x1-x)./(x1-x0);
    aux= zeros(2,2);
    for j = 1:2
        for k =1:2
            aux(j,k) = (x1-x0) * dot(fy{j}(x0+(x1-x0)*node).*fy{k}(x0+(x1-x0)*node),W);
        end
    end
    Mass(t(i,:),t(i,:)) = Mass(t(i,:),t(i,:))+aux;
end
Mass =Mass(nf,nf);