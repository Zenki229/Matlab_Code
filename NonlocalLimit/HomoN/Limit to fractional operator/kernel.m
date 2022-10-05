        
function [ v ] = kernel(s, p, del)
        c = 4^p*p*gamma(p+1/2)/(pi^(1/2)*gamma(1-p));
        v = c*s.^(-1-2*p).*(s<=del & s>=0);
%         v = exp(-abs(s)).*(s<delta & s>0) / (1-exp(-delta));
end