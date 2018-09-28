%This function describes the motion of quadrotor as translation motion in
%x y z axis

clc; clear;

Time = 20;
dt = 0.001;
t = 0:dt:Time;
N = length(t);

% constants
m = 2; 
g = 9.81;
L = 0.08;
I = diag([0.081, 0.0812, 0.1320]);
Jr = 5 * 10^-5;

% of the blade
area = 10^-3; 
radius = 0.03;

%Aerodynamic constants
rho = 1;
Cd = 0.19; Ct = 0.23;

Kf = 0.5 * rho * area * Ct * radius * radius;
Km = 0.5 * rho * area * Cd * radius * radius;
Kd = diag([0.7 0.7 1.4]) * 10^-4;
matrix = [Kf Kf Kf Kf;0 Kf*L 0 -Kf*L;-Kf*L 0 Kf*L 0;Km -Km Km -Km];

% state
x = zeros(3,N);
v = zeros(3,N);
R = zeros(3,3,N);
phi = zeros(1,N);
the = zeros(1,N);
psi = zeros(1,N);
omega = zeros(3,N);
w = [6884.13; 2000; 6884.13; 0];

wp = w ./ 100;

% initialize
x(:,1) =[5; 5; 5];
v(:,1)=[0; 0; 0];
R(:,:,1) = vrrotvec2mat([0,0,1,90*pi/180]);
omega(:,1) = zeros(3,1);
omega(:,1) = [0;0;0];

T = m*g;
M = zeros(3,N);
e3=[0; 0; 1];

I_inv = inv(I);

% Controller variables
phi_des = zeros(1,N);
the_des = zeros(1,N);
psi_des = zeros(1,N);

phi_des(1) = 0;
the_des(1) = 0;
psi_des(1) = 0;
p_des = 0;
q_des = 0;
r_des = 0;

kp_pos = 10;
kd_pos = 12;

kp_att = 75;
kd_att = 5;

x_des = 0;  % Reach from 0 in Time sec
y_des = 0;
z_des = 0;

vx = 0;
vy = 0;
vz = 0;

for n=1:N-1

    
    u1 = m * [0; 0; 0] + kd_pos * [vx - v(1,n); vy - v(2,n); vz - v(3,n)]...
        + kp_pos * [x_des - x(1,n); y_des - x(2,n); z_des - x(3,n)] + m * g * e3;
    
    T(n) = e3' * R(:,:,n)' * u1;
    
    bz_des = u1 / norm(u1);
    bz = R(:,3,n);
    angle_err = inv(R(:,:,n)) * Hat(bz) * bz_des;
    
    phi_des(n+1) = angle_err(1);
    the_des(n+1) = angle_err(2);
    psi_des(n+1) = 0;
    
    p_des = (phi_des(n+1) - phi_des(n)) / dt;
    q_des = (the_des(n+1) - the_des(n)) / dt;
    r_des = (psi_des(n+1) - psi_des(n)) / dt;
    
    [phi(n), the(n), psi(n)] = GetEulerAngles(R(:,:,n));
    p = omega(1,n);
    q = omega(2,n);
    r = omega(3,n);
    
    M(1,n) = kp_att * (phi_des(n) - phi(n)) + kd_att * (p_des - p);
    M(2,n) = kp_att * (the_des(n) - the(n)) + kd_att * (q_des - q);
    M(3,n) = kp_att * (psi_des(n) - psi(n)) + kd_att * (r_des - r);
    
    omega(:,n+1) = omega(:,n) + dt * I_inv * (M(:,n) - Hat(omega(:,n)) * I * omega(:,n));
    R(:,:,n+1) = R(:,:,n) * expm(Hat(omega(:,n)) * dt);
    v(:,n+1) = v(:,n) + dt*(-m * g * e3 + (T(n) * R(:,:,n) * e3)) / m;
    x(:,n+1) = x(:,n) + v(:,n) * dt;
end

%% Animation using hgtransform
% 
% view(3);
% grid on;
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% set(gca,'XLim',[-5 15],'YLim',[-5 15],'ZLim',[-10 10]);
% 
% [x1 y1 z1] = cylinder([0.2 0.2]);
% [x2 y2 z2] = cylinder([0.2 0.2]);
% g(1) = surface((z1-0.5)*2,y1,x1,'FaceColor','red');
% g(2) = surface(x2,(z2-0.5)*2,y2,'FaceColor','blue');
% 
% [h1 y3 z3] = cylinder([0.1 0.1]);
% [x4 y4 z4] = cylinder([0.1 0.1]);
% a(1) = surface(h1-1,y3*2+0.3,z3*0.2+0.1,'FaceColor','red');
% a(2) = surface(x4-1,y4*2-0.3,z4*0.2+0.1,'FaceColor','red');
% a(3) = surface(h1-1,y3,z3*0.2+0.1,'FaceColor','yellow');
% 
% c(1) = surface(h1+1,y3*2+0.3,z3*0.2+0.1,'FaceColor','red');
% c(2) = surface(x4+1,y4*2-0.3,z4*0.2+0.1,'FaceColor','red');
% c(3) = surface(h1+1,y3,z3*0.2+0.1,'FaceColor','yellow');
% 
% b(1) = surface(h1,y3*2+1.3,z3*0.2+0.1,'FaceColor','blue');
% b(2) = surface(x4,y4*2+0.7,z4*0.2+0.1,'FaceColor','blue');
% b(3) = surface(h1,y3+1,z3*0.2+0.1,'FaceColor','yellow');
% 
% d(1) = surface(h1,y3*2-0.7,z3*0.2+0.1,'FaceColor','blue');
% d(2) = surface(x4,y4*2-1.3,z4*0.2+0.1,'FaceColor','blue');
% d(3) = surface(h1,y3-1,z3*0.2+0.1,'FaceColor','yellow');
% 
% quad = hgtransform;
% set(g, 'parent', quad);
% 
% r1 = hgtransform;
% set(a, 'parent', r1);
% 
% r2 = hgtransform;
% set(b, 'parent', r2);
% 
% r3 = hgtransform;
% set(c, 'parent', r3);
% 
% r4 = hgtransform;
% set(d, 'parent', r4);
% drawnow
% 
% T1x = makehgtform('translate',[-1 0 0]);
% T2x = makehgtform('translate',[1 0 0]);
% 
% T1y = makehgtform('translate',[0 -1 0]);
% T2y = makehgtform('translate',[0 1 0]);
% 
% latitude = x(1,:);
% longitude = x(2,:);
% altitude = x(3,:);
% [axis,theta] = axang(R);
% ang = [0; 0; 0; 0];
% 
% for i = 1:length(latitude)
%     translation = makehgtform('translate',latitude(i),longitude(i),altitude(i));
%     rotation = makehgtform('axisrotate',axis(:,i),theta(i));
%     
%     zb = R(:,:,i) * [0;0;1];
%     R1 = makehgtform('axisrotate',zb,ang(1));
%     R1(1,4) = latitude(i);
%     R1(2,4) = longitude(i);
%     R1(3,4) = altitude(i);
%     
%     R2 = makehgtform('axisrotate',zb,ang(2));
%     R2(1,4) = latitude(i);
%     R2(2,4) = longitude(i);
%     R2(3,4) = altitude(i);
%     
%     R3 = makehgtform('axisrotate',zb,ang(3));
%     R3(1,4) = latitude(i);
%     R3(2,4) = longitude(i);
%     R3(3,4) = altitude(i);
%     
%     R4 = makehgtform('axisrotate',zb,ang(4));
%     R4(1,4) = latitude(i);
%     R4(2,4) = longitude(i);
%     R4(3,4) = altitude(i);
%     
%     c1 = R(:,:,i) * [-1; 0; 0];
%     c2 = R(:,:,i) * [0; 1; 0];
%     c3 = R(:,:,i) * [1; 0; 0];
%     c4 = R(:,:,i) * [0; -1; 0];
%     
%     x1 = makehgtform('translate', [1 0 0]);
%     x2 = makehgtform('translate', [-1 0 0]);
%     y1 = makehgtform('translate', [0 1 0]);
%     y2 = makehgtform('translate', [0 -1 0]);
%     
%     h1 = makehgtform('translate', [c1(1,1) c1(2,1) c1(3,1)]);
%     h2 = makehgtform('translate', [c2(1,1) c2(2,1) c2(3,1)]);
%     h3 = makehgtform('translate', [c3(1,1) c3(2,1) c3(3,1)]);
%     h4 = makehgtform('translate', [c4(1,1) c4(2,1) c4(3,1)]);
%     
%     set(quad, 'matrix', translation*rotation);
%     set(r1, 'matrix', h1*R1*rotation*x1);
%     set(r2, 'matrix', h2*R2*rotation*y2);
%     set(r3, 'matrix', h3*R3*rotation*x2);
%     set(r4, 'matrix', h4*R4*rotation*y1);
%    
%     ang = ang + (wp .* dt);
%     
%     p0 = [latitude(i);longitude(i);altitude(i)];
%     p1 = p0 + zb;
%     hold on;
%     plot3([p0(1), p1(1)], [p0(2), p1(2)], [p0(3), p1(3)]);
%     plot3(latitude,longitude,altitude);
%     hold on;
%     pause(0.01);
% end

%% plotting the 3D graph as an animation

% curve = animatedline('LineStyle',':','Marker','o','MarkerFaceColor','r');
% set(gca,'XLim',[-4.5 0],'YLim',[-1 1],'ZLim',[-2 2]);
% view(43,24);
% grid on;
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% for i = 1:length(t)
%     clearpoints(curve);
%     addpoints(curve,x(1,i),x(2,i),x(3,i));
%     drawnow
% end

%% Ploting normal graphs

figure

s = t(1,1:20000);

subplot(2,2,1);
plot(t,x(1,:));
grid on;
title('x')
xlabel('time (s)'); ylabel('x (m)');

hold on;
plot(s,x_des,'r--');


subplot(2,2,2);
plot(t,x(2,:));
grid on;
title('y')
xlabel('time (s)'); ylabel('y (m)');

hold on;
plot(s,y_des,'r--');

subplot(2,2,3);
plot(t,x(3,:));
grid on;
title('z')
xlabel('time (s)'); ylabel('z (m)');

hold on;
plot(s,z_des,'r--');

subplot(2,2,4);
plot(x(1,:),x(2,:));
grid on;
title('xy')
xlabel('x (m)'); ylabel('y (m)');

hold on;
plot(x_des,y_des,'r--');

figure
subplot(2,2,1);
plot(t,phi);
grid on;
title('phi')
xlabel('time (s)'); ylabel('phi (rad)');

hold on;
plot(t,phi_des,'r--');

subplot(2,2,2);
plot(t,the);
grid on;
title('theta')
xlabel('time (s)'); ylabel('theta (rad)');

hold on;
plot(t,the_des,'r--');

subplot(2,2,3);
plot(t,psi);
grid on;
title('psi')
xlabel('time (s)'); ylabel('psi (rad)');

hold on;
plot(t,psi_des,'r--');

figure
subplot(2,2,1);
plot(s,T);
grid on;
title('Thrust')
xlabel('time (s)'); ylabel('T (N)');

subplot(2,2,2);
plot(t(2:end),M(1,2:end));
grid on;
title('Mx')
xlabel('time (s)'); ylabel('Mx (Nm)');

subplot(2,2,3);
plot(t(2:end),M(2,2:end));
grid on;
title('My')
xlabel('time (s)'); ylabel('My (Nm)');

subplot(2,2,4);
plot(t(2:end),M(3,2:end));
grid on;
title('Mz')
xlabel('time (s)'); ylabel('Mz (Nm)');


