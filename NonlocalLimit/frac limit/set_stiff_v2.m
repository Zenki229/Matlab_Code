function [ c ] = set_stiff_v2( al, M, h )
%SET_STIFF Summary of this function goes here
%   Detailed explanation goes here

M
h
c=zeros(M-1,1);

c(1)=local_stiff_1(al,h);
c(2)=local_stiff_2(al,h);
for i=3:M-1
    c(i)=local_stiff_3(i,al,h);    
end
%Stiff=Stiff+Stiff'-diag(diag(Stiff));


end



