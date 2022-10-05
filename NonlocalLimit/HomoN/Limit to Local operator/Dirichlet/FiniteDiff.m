function [A,h,xx, u ,error] = FiniteDiff( h, r, f ,g)
%%%%%%%% input %%%%%%%%%%
%%%  h:mesh size  %%%
%%%  r=horizon/h  %%%
%%%  f: right hand side %%%
%%%  g: boundary condition (just input exact solution u) %%%

%%%%%%%% output %%%%%%%%%%%
%%%  A: stiffness matrix %%%%
%%%  u: numerical solution at grid points %%%
%%%  error: max error at grid points %%%

%%%%%%%% example usage %%%%%%
%%%  [A, u ,error] = FiniteDiff( 2^-5, 3, @f2 ,@u2) %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N=1/h-1; xx=0:h:1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
A=eye(N)*4*quad(@(x)ga(r*h,x),0,r*h,10^-12); %method 1
for i=1:r
    B=[zeros(i,N);-2*quad(@(x)ga(r*h,x),(i-1)*h,i*h,10^-12)*eye(N-i,N)];
    A=A+B+B';
end
t=h:h:N*h;
F=f(t');

%boundary condition
for i=1:r
    s1=0;s2=0;
    for j=r:-1:i
        t=quad(@(x)ga(r*h,x),(j-1)*h,j*h,10^-12);
        s1=s1+g((i-j)*h)*2*t;
        s2=s2+g((N-i+j+1)*h)*2*t;
    end
     F(i)=F(i)+s1;
     F(N-i+1)=F(N-i+1)+s2;
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
A=zeros(N); s=0; %method 2
for i=1:r  
    t=quad(@(x)ga(r*h,x).*x,(i-1)*h,i*h,10^-14)/(i*h);
    s=s+4*t;
    B=[zeros(i,N);-2*t*eye(N-i,N)];
    A=A+B+B';
end
A=A+eye(N)*s;
t=h:h:N*h;
F=f(t');

%boundary condition
for i=1:r
    s1=0;s2=0;
    for j=r:-1:i
        t=quad(@(x)ga(r*h,x).*x,(j-1)*h,j*h,10^-14)/(j*h);
        s1=s1+g((i-j)*h)*2*t;
        s2=s2+g((N-i+j+1)*h)*2*t;
    end
     F(i)=F(i)+s1;
     F(N-i+1)=F(N-i+1)+s2;
end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%{
A=zeros(N); s=0; %method 3 (convergent method for fixed r)
for i=1:r  
    t=quad(@(x)ga(r*h,x).*x.^2,(i-1)*h,i*h,10^-14)/(i*h)^2;
    s=s+4*t;
    B=[zeros(i,N);-2*t*eye(N-i,N)];
    A=A+B+B';
end   
A=A+eye(N)*s;
t=h:h:N*h;
F=f(t');

%boundary condition
for i=1:r
    s1=0;s2=0;
    for j=r:-1:i
        t=quad(@(x)ga(r*h,x).*x.^2,(j-1)*h,j*h,10^-14)/(j*h)^2;
        s1=s1+g((i-j)*h)*2*t;
        s2=s2+g((N-i+j+1)*h)*2*t;
    end
     F(i)=F(i)+s1;
     F(N-i+1)=F(N-i+1)+s2;
end
 
%%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u=A\F;

%numerical error
error=0;
for i=1:N
    error=max(error, abs(u(i)-g(i*h)));
end
end