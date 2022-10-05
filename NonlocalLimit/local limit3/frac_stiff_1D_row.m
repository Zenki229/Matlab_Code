function [c] = frac_stiff_1D_row(h,p,delta,kerf)
        N = 1/h - 1;
        r = floor(delta / h+0.0001);
% ignore r<=1

 %m = |i-j| = 0 free element
        ann1 = quadgk(@(t)kerf(t, p, delta).*(  2*t.^2 / h^2), 0, h);
        ann2 = quadgk(@(t)kerf(t, p, delta).*( - 2*t.^2 / h^2 +...
            4*t / h - 4/3), h, 2*h);
        ann3 = quadgk(@(t)kerf(t, p, delta), 2*h, delta);
        ele_free = ann1 + ann2 + 4*ann3 / 3;
 %m = |i-j| = 0 bdry element     
        ann1 = quadgk(@(t)kerf(t, p, delta).*( 2*t.^2 / h^2), 0, h);
        ann2 = 0;
        ann3 = quadgk(@(t)kerf(t, p, delta), h, delta);
        ele_bdry = ann1/2 + ann2 + 2*ann3 / 3;
        %m = |i-j| = 1
        anplus1 = quadgk(@(t)kerf(t, p, delta).*( - t.^2 / h^2), 0, h);
        anplus2 = quadgk(@(t)kerf(t, p, delta).*(...
            5*t.^2 / (2*h^2) - 7*t / (2*h) + 7 / 6), h, 2*h);
        anplus3 = quadgk(@(t)kerf(t, p, delta).*( - 3*t.^2 / (2*h^2) +...
            9*t / (2*h) - 25 / 6), 2*h, 3*h);
        anplus4 = quadgk(@(t)kerf(t, p, delta), 3*h, delta) / 3;
        anplus = anplus1 + anplus2 + anplus3 + anplus4;
        
        % set r>=2
        andpls = zeros(r-1, 1);
        for m=2:r+1
            andpls1 = quadgk(@(t)kerf(t, p, delta).*( ...
                (m-2)*t.^2 / (2*h^2) - (m^2-4*m+4)*t / (2*h) +...
                (m^3-6*m^2+12*m-8) / 6), (m-2)*h, (m-1)*h);
            andpls2 = quadgk(@(t)kerf(t, p, delta).*( -...
                (3*m-2)*t.^2 / (2*h^2) + (3*m^2-4*m)*t / (2*h) -...
                (3*m^3-6*m^2+4) / 6), (m-1)*h, (m)*h); 
            andpls3 = quadgk(@(t)kerf(t, p, delta).*(...
                (3*m+2)*t.^2 / (2*h^2) - (3*m^2+4*m)*t / (2*h) +...
                (3*m^3+6*m^2-4) / 6), (m)*h, (m+1)*h);
            andpls4 = quadgk(@(t)kerf(t, p, delta).*( -...
                (m+2)*t.^2 / (2*h^2) + (m^2+4*m+4)*t / (2*h) -...
                (m^3+6*m^2+12*m+8) / 6), (m+1)*h, (m+2)*h);
            andpls(m-1) = andpls1 + andpls2 + andpls3 + andpls4;
        end
        

        for kj = 0:r+1
            if kj == 0 
               c(kj+1)=ele_free;
            elseif kj == 1
                c(kj+1)=anplus;
            else
                c(kj+1)=andpls(kj-1);
            end
        end    
        c=h*[c,zeros(1,N-length(c))];

    end