function S = Stiff_nonlocal_mat_v2(row,h,delta)
global xb1 xb2 x_int;
    r = floor(delta/h+1e-5);
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
% Contracti