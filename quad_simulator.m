%This function describes the motion of quadrotor as translation motion in
%x y z axis

clc; clear;

Time = 5;
dt = 0.001;
t = 0:dt:Time;
N = length(t);

% constants
m = 2; g = 9.81;
kf=3.2*10^-6*.12;
L=.5;
I=diag([0.081, 0.0812, 0.1320]);
b=3.2*10^-6;

% state
x = zeros(3,N);
v = zeros(3,N);
R = zeros(3,3,N);
omega = zeros(3,N);

% initialize
x(:,1) =[0; 0; 0];
v(:,1)=[0; 0; 0];
R(:,:,1) = vrrotvec2mat([1,0,0,2*pi/180]);
omega(:,1) = zeros(3,1);

T=m*g;
M = zeros(3,N);
e3=[0; 0; 1];

I_inv = inv(I);

for n=1:N-1
    omega(:,n+1) = omega(:,n) + dt*I_inv*(M(:,1)-Hat(omega(:,n))*I*omega(:,n));
    R(:,:,n+1) = R(:,:,n)*expm(Hat(omega(:,n))*dt);
    v(:,n+1) = v(:,n) + dt*(m*g*e3 - (T*R(:,:,n)*e3))/m;
    x(:,n+1) = x(:,n) + v(:,n)*dt;
end

%% plot

figure
plot(t,x);
grid on;
title('x')
xlabel('time (s)'); ylabel('x (m)');

figure
plot(t,v);
grid on;
title('v')
xlabel('time (s)'); ylabel('v (m/s)');

figure
plot(t,omega);
grid on;
title('\omega')
xlabel('time (s)'); ylabel('\omega (rad/s)');
