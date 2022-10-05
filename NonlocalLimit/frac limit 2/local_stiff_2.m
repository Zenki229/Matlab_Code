function [ SL ] = local_stiff_2(al,h)
%LOCAL_STIFF_2 Summary of this function goes here
%   Detailed explanation goes here

a22=-2*h^(1-al)/(2-al)/(3-al);
%a12=-quad2d(@(x,y) (x+y).*y.*(y-x).^(-al-1),-h,0,0,h)/(h^2);
a12=-integral2(@(x,y) (x+y).*y.*(y-x).^(-al-1),-h,0,0,h)/(h^2);
a23=a12;
%a13=quad2d(@(x,y)x.*y.*(3*h+y-x).^(-al-1),0,h,-h,0)/(h^2);
a13=integral2(@(x,y)x.*y.*(3*h+y-x).^(-al-1),0,h,-h,0)/(h^2);
aout=quadgk(@(x) (2*h-x).*(x-h).*(x).^(-al),h,2*h)/(h^2)*2/al;
SL=a22+2*(a12+a23+a13+aout);

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


