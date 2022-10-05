        
function [ v ] = kernel(s, p, del)
        v = (2-2*p)*s.^(-1-2*p).*(s<=del & s>=0)* del^(2*p-2)/2;
%         v = exp(-abs(s)).*(s<delta & s>0) / (1-exp(-delta));
end