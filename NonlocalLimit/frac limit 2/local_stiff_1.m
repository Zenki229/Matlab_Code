function [ SL ] = local_stiff_1(al,h)
%LOCAL_STIFF_1 Summary of this function goes here
%   Detailed explanation goes here


%a=( quadgk(@(x) x.^(2-al)+x.^2.*(2*h-x).^(-al),0,h) + quadgk(@(x) x.^(-al).*(2*h-x).^2+(2*h-x).^(2-al), h, 2*h))/al/h^2;
a=4*quadgk(@(x) x.^(2-al)+x.^2.*(2*h-x).^(-al),0,h)/al/h^2;
%b= 4*h^(1-al)/(2-al)/(3-al) + 2*quad2d(@(x,y) (x+y).^2.*(y-x).^(-al-1),-h,0,0,h)/(h^2);
b= 4*h^(1-al)/(2-al)/(3-al) + 2*integral2(@(x,y) (x+y).^2.*(y-x).^(-al-1),-h,0,0,h)/(h^2);

SL=a+b;

end

