function [ vh ] = proj_l2(v,xx)

N=length(xx)-2;
vh=zeros(N,1);

for j=1:N
    
    x0=xx(j);
    x1=xx(j+1);
    x2=xx(j+2);
    h1=x1-x0;
    h2=x2-x1;
    phi=@(x) 0+(x-x0)./h1.*(x>=x0 & x<=x1)+(1-(x-x1)./h2).*(x>x1 & x<=x2);
    FF=@(x) v(x).*phi(x);
    
    vh(j)= quadgk(FF,x0,x2);
end


%vh=M\vh;

end

