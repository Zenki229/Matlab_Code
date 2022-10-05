function [ c ] = Stiff_nonlocal_row_free(h, p, delta, kerf)
        N = 1/h - 1;
        r = floor(delta / h+0.0001);

        
        if r<=1
            an10 =  quadgk(@(t)kerf(t, p, delta).*t.^2, 0, h);
            an1 = 2*an10 / h^2;
            an20 =  quadgk(@(t)kerf(t, p, delta).*t.^3, 0, h);
            an2 = an20 / h^3;
            an = an1 - an2;
            an30 =  quadgk(@(t)kerf(t, p, delta).*t.^3, 0, h);
            an3 = 2*an30 / (3*h^3);
            an40 =  quadgk(@(t)kerf(t, p, delta).*t.^2, 0, h);
            an4 = an40 / h^2;      
            anpls = an3 - an4;
            an50 =  quadgk(@(t)kerf(t, p, delta).*t.^3, 0, h);
            an5 = -an50 / (6*h^3);  
%             A1 = diag(an*ones(N, 1));
%             A2 = diag(anpls*ones(N-1, 1), 1);
%             A3 = diag(an5*ones(N-2, 1), 2);
%             A = A1 + A2 + A2' + A3 + A3';
%             SN = sparse(A);
            c=h*[an,anpls,an5,zeros(1,N-3)];
            return;
        end
        
        %m = |i-j| = 0 
        ann1 = quadgk(@(t)kerf(t, p, delta).*(-t.^3 / h^3 + 2*t.^2 / h^2), 0, h);
        ann2 = quadgk(@(t)kerf(t, p, delta).*(t.^3 / (3*h^3) - 2*t.^2 / h^2 +...
            4*t / h - 4/3), h, 2*h);
        ann3 = quadgk(@(t)kerf(t, p, delta), 2*h, delta);
        ann = ann1 + ann2 + 4*ann3 / 3;
        
        %m = |i-j| = 1
        anplus1 = quadgk(@(t)kerf(t, p, delta).*(2*t.^3 / (3*h^3) - t.^2 / h^2), 0, h);
        anplus2 = quadgk(@(t)kerf(t, p, delta).*(-t.^3 / (2*h^3) +...
            5*t.^2 / (2*h^2) - 7*t / (2*h) + 7 / 6), h, 2*h);
        anplus3 = quadgk(@(t)kerf(t, p, delta).*(t.^3 / (6*h^3) - 3*t.^2 / (2*h^2) +...
            9*t / (2*h) - 25 / 6), 2*h, 3*h);
        anplus4 = quadgk(@(t)kerf(t, p, delta), 3*h, delta) / 3;
        anplus = anplus1 + anplus2 + anplus3 + anplus4;
        
        % set r>=2
        andpls = zeros(r-1, 1);
        for m=2:r+1
            andpls1 = quadgk(@(t)kerf(t, p, delta).*(-t.^3 / (6*h^3) +...
                (m-2)*t.^2 / (2*h^2) - (m^2-4*m+4)*t / (2*h) +...
                (m^3-6*m^2+12*m-8) / 6), (m-2)*h, (m-1)*h);
            andpls2 = quadgk(@(t)kerf(t, p, delta).*(t.^3 / (2*h^3) -...
                (3*m-2)*t.^2 / (2*h^2) + (3*m^2-4*m)*t / (2*h) -...
                (3*m^3-6*m^2+4) / 6), (m-1)*h, (m)*h); 
            andpls3 = quadgk(@(t)kerf(t, p, delta).*(-t.^3 / (2*h^3) +...
                (3*m+2)*t.^2 / (2*h^2) - (3*m^2+4*m)*t / (2*h) +...
                (3*m^3+6*m^2-4) / 6), (m)*h, (m+1)*h);
            andpls4 = quadgk(@(t)kerf(t, p, delta).*(t.^3 / (6*h^3) -...
                (m+2)*t.^2 / (2*h^2) + (m^2+4*m+4)*t / (2*h) -...
                (m^3+6*m^2+12*m+8) / 6), (m+1)*h, (m+2)*h);
            andpls(m-1) = andpls1 + andpls2 + andpls3 + andpls4;
        end
        

        for kj = 0:r+1
            if kj == 0
%                 Ak = diag(ann*ones(N,1));
               c(kj+1)=ann;
            elseif kj == 1
                %Ak = diag(anplus*ones(N-1, 1), 1);
                c(kj+1)=anplus;
%                 Ak = Ak + Ak';
            else
                %Ak = diag(andpls(kj-1)*ones(N-kj, 1), kj);
                c(kj+1)=andpls(kj-1);
%                 Ak = Ak + Ak';
            end
%             SN = SN + Ak;
        end    
        %SN = sparse(SN);
        c=h*[c,zeros(1,N-length(c))];

    end