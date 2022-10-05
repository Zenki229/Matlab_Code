function [c] = Stiff_nonlocal_collocation(del,h,s,kerfun)

r = round(del/h+0.000001);
if r<=1
    a1 = -1/h^2 * quadgk(@(x) x.^2.*kerfun(x,s,del),0,del);
    c = [-2*a1;a1];
    return
end
cc = zeros(r,1);
for m = 1:r
    cc(m) =  -1/(h*m)^2* quadgk(@(x) x.^2.*kerfun(x,s,del),(m-1)*h,m*h);
end
aux = sum(cc)*(-2);
c = [aux;cc];