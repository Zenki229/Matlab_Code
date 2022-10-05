function S = Stiff_nonlocal_mat_2(row,h,delta)
    r = delta/h;
    xb1 = -delta+h:h:0;
    xb2 = 1:h:1+delta-h;
    x_int = h:h:1-h;
    nn1 = numel(xb1);
    nn2 = numel(xb2);
    nn = numel(x_int);
    nn_all = nn1+nn+nn2;
    S = zeros(nn_all,nn_all);
for i = 1:nn_all
        for j = 1:numel(row)
            if i+j-1>nn_all 
                S(i,i) = S(i,i)+row(j);
            else 
            S(i,i+j-1) = S(i,i+j-1)+row(j);
            end
        end
    end
    for i = 1:nn_all
        for j = 2:numel(row)
            if i-j+1 < 1 
                S(i,i) = S(i,i)+row(j);
            else 
            S(i,i-j+1) = S(i,i-j+1)+row(j);
            end
        end
    end