function [ c ] = set_stiff( M, al )
%SET_STIFF Summary of this function goes here
%   Detailed explanation goes here
h=1/M;
c=zeros(M-1,1);

M

h

c(1)=local_stiff_1(al,h);
c(2)=local_stiff_2(al,h);
for i=3:M-1
    c(i)=local_stiff_3(i,al,h);    
end
%Stiff=Stiff+Stiff'-diag(diag(Stiff));


end



