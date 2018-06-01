%This function describes the motion of quadrotor as translation motion in
%x y z axis
function [x,t,v]=integrate(N,h,m)
v(1,:)=[10 10 0];
x(1,:) =[5 6 0];
g=9.81;
T=m*g;
e3=[0 0 1];
R=diag([1 1 1]);
for n=1:N
    v(n+1,:)=v(n,:)+h*(m*g*e3-T*e3*R)/m;
    x(n+1,:)=x(n,:)+v(n,:)*h;
   t(n+1)=n*h;
    %plot(v(n+1),x(n+1));
    %hold on;
end
%R(n+1)=R(n)*expm((wn-wn')*h);
%w(n+1)=w(n)+h*

   