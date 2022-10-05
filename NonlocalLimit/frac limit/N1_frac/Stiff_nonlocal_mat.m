function S = Stiff_nonlocal_mat(row,h,delta)
    r = delta/h;
    xb1 = -delta+h:h:0;
    xb2 = 1:h:1+delta-h;
    x_int = h:h:1-h;
    nn1 = numel(xb1);
    nn2 = numel(xb2);
    nn = numel(x_int);
    nn_all = nn1+nn+nn2;
    S = zeros(nn_all,nn_all);
    for i = 1:nn1
        for j = 1:nn
                if nn1-i+1+j <= numel(row)
                    S(i,nn1+j) =  row(nn1-i+1+j);
                end
        end
        S(i,i) = -sum(S(i,:));
    end
    for i = nn1+(1:nn)
        for j = 1:numel(row)
            if i+j-1>nn_all 
                S(i,i) = S(i,i)+row(j);
            else 
            S(i,i+j-1) = S(i,i+j-1)+row(j);
            end
        end
    end
    for i = nn1+(1:nn)
        for j = 2:numel(row)
            if i-j+1 < 1 
                S(i,i) = S(i,i)+row(j);
            else 
            S(i,i-j+1) = S(i,i-j+1)+row(j);
            end
        end
    end
    for i = nn_all-nn2+1:nn_all
        for j = 1:nn
                if nn1-(nn_all-i+1)+1+j <= numel(row)
                    S(i,nn1+nn+1-j) = row(nn1-(nn_all-i+1)+1+j);
                end
        end
        S(i,i) = -sum(S(i,:));
    end