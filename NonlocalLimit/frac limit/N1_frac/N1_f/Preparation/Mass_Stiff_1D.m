function [Mass,Stiff]=Mass_Stiff_1D(vert,t,nf)
% generate mass and stiff mat with nonzero bdry condition

nn = size(vert,1);
nt = size(t,1);
Mass = zeros(nn,nn);
Stiff = zeros(nn,nn);
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
%             aux(j,k) = quadgk(@(x) fy{j}(x).*fy{k}(x),x0,x1);
            aux(j,k) = (x1-x0) * dot(fy{j}(x0+(x1-x0)*node).*fy{k}(x0+(x1-x0)*node),W);
        end
    end
%     dfy = [1/(x1-x0);-1/(x1-x0)];
    bux = 1/(x1-x0)*[1 -1;-1 1];
    Mass(t(i,:),t(i,:)) = Mass(t(i,:),t(i,:))+aux;
    Stiff(t(i,:),t(i,:)) = Stiff(t(i,:),t(i,:))+bux;
end
Mass =Mass(nf,nf);
Stiff=Stiff(nf,nf);

