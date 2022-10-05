function [Mass] = Weighted_Mass_1D(vert,t,nf,q)
% generate mass and stiff mat with nonzero bdry condition

nn = size(vert,1);
nt = size(t,1);
Mass = zeros(nn,nn);
[node,W] = gauleg(0,1,20);
for i = 1:nt
    x0 = vert(t(i,1));
    x1 = vert(t(i,2));
    fy{1} = @(x) (x-x0)./(x1-x0);
    fy{2} = @(x) (x1-x)./(x1-x0);
    aux= zeros(2,2);
    for j = 1:2
        for k =1:2
            u = @(x) q((x1-x0).*x+x0).*fy{j}((x1-x0).*x+x0).*fy{k}((x1-x0).*x+x0);
            aux(j,k) = (x1-x0)*dot(u(node),W);
        end
    end
    Mass(t(i,:),t(i,:)) = Mass(t(i,:),t(i,:))+aux;

end
Mass =Mass(nf,nf);




