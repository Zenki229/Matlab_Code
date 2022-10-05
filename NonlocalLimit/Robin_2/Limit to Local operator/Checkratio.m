function [ Ratio ] = Checkratio(v, p)

n=length(v)-1;

for j=1:n
    Ratio(j)=v(j)/v(j+1);
end


Ratio=log(Ratio)./log(p(2:end)./p(1:end-1));



end