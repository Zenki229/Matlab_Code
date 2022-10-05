function [ SL ] = local_stiff_3(k,al,h)
%LOCAL_STIFF_2 Summary of this function gokes here
%   Detailed explanation goes here



a13=-integral2(@(x,y) x.*y.*((k-1)*h+y-x).^(-al-1),0,h,0,h)/(h^2);
a14=integral2(@(x,y) x.*y.*((k+1)*h+y-x).^(-al-1),0,h,-h,0)/(h^2);
a23=integral2(@(x,y) x.*y.*((k-3)*h+y-x).^(-al-1),-h,0,0,h)/(h^2);
%a24=-quad2d(@(x,y) x.*y.*((k-1)*h+y-x).^(-al-1),-h,0,-h,0)/(h^2);
a24=a13;

SL=2*(a13+a14+a23+a24);
% 
% if i==j
%     SL=1/al*(quadgk(@(x) phi((x-x(i))/h).^2.*(x.^(-al)+(1-x).^(-al)), x(i),x(i+1))+...
%         + quadgk(@(x) psi((x(i+2)-x)/h).^2.*(x.^(-al)+(1-x).^(-al)), x(i+1),x(i+2)));
% else if j-i==1
%         SL=;
%     else
%         SL=;      
%     end
% end
end