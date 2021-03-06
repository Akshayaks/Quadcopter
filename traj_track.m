%This function describes the motion of quadrotor as translation motion in
%x y z axis

clc; clear;

Time = 15;
dt = 0.01;
t = 0:dt:Time;
N = length(t);

% constants
m = 2; 
g = 9.81;
L = 0.08;
I = diag([0.081, 0.0812, 0.1320]);
Jr = 5 * 10^-5;

% %of the blade
area = 10^-3; 
radius = 0.03;
L = 0.08; %distance from the axis of rotation to the rotors

%Aerodynamic constants
rho = 1;
Cd = 0.19; Ct = 0.23;

% Kp and Kd for phi theta psi in this order
% K = [141.0, 20.0; 5.5, 1.5; 141.0, 20.0];

K = [250.0, 20.0; 4, 0.5; 141.0, 20.0];


% F = Kf*omega*omega;
% M = Km*omega*omega;
% M = matrix * omega;
% Kf = 3.2*10^-6 * 0.12;
% Km = 0.7 * 10^-4;

Kf = 0.5 * rho * area * Ct * radius * radius;
Km = 0.5 * rho * area * Cd * radius * radius;
Kd = diag([0.7 0.7 1.4]) * 10^-4;
matrix = [Kf Kf Kf Kf;0 Kf*L 0 -Kf*L;-Kf*L 0 Kf*L 0;Km -Km Km -Km];

% state
x = zeros(3,N);
v = zeros(3,N);
phi = zeros(1,N);
the = zeros(1,N);
psi = zeros(1,N);
p = 0;
q = 0;
r = 0;
p_dot = 0;
q_dot = 0;
r_dot = 0;

R = zeros(3,3,N);
omega = zeros(3,N);
w = [6195.717; 6884.13; 6884.13; 6884.13];

wp = w ./ 100;

% initialize
x(:,1) =[0; 0; 5];
v(:,1)=[0; 0; 0];
R(:,:,1) = vrrotvec2mat([0,0,1,90*pi/180]);
phi(1) = 0;
the(1) = 0;
psi(1) = 1.5708;
omega(:,1) = zeros(3,1);
omega_dot = [0; 0; 0];


T = m*g;
M = zeros(3,N);
e3=[0; 0; 1];

I_inv = inv(I);

% Getting the trajectory along x, y and z directions: Gives all desired
% values as a function of time

rx = traj_gen(0,10,Time)';  % Reach from a to b in 5 sec
ry = traj_gen(0,10,Time)';
rz = 5;

vx = polyder(rx);
vy = polyder(ry);
vz = 0;

ax = polyder(vx);
ay = polyder(vy);
az = 0; 

phi_des = zeros(1,N);
theta_des = zeros(1,N);
psi_des = zeros(1,N);

psi_des(1) = 1.5708;
p_des = 0;
q_des = 0;
r_des = 0;

for n=1:N-1
    [phi(n+1), the(n+1), psi(n+1)] = GetEulerAngles(R(:,:,n));
    phi_des(n+1) = (1/g) * (polyval(ax,n*dt) * sin(psi(n)) - polyval(ay,n*dt) * cos(psi(n)));
    theta_des(n+1) = (1/g) * (polyval(ax,n*dt) * cos(psi(n)) + polyval(ay,n*dt) * sin(psi(n)));
    psi_des(n+1) = 1.5708;
    if n == 1
        ;
    else
        p = (phi(n+1) - phi(n)) / dt;
        q = (the(n+1) - the(n)) / dt;
        r = (psi(n+1) - psi(n)) / dt;
        p_des = (phi_des(n+1) - phi_des(n))/dt; 
        q_des = (theta_des(n+1) - theta_des(n))/dt;
        r_des = (psi_des(n+1) - psi_des(n))/dt;
    end
    
    
    T_des = m * g + m * az;
    M_des = I * omega_dot + Hat(omega(:,n)) * I * omega(:,n);
    
    w = sqrt(inv(matrix) * [T_des;M_des]);
    
    % Theta double dot = kp*(error in theta) + kv*(error in theta dot)
    omega_dot = [K(1,1) * (phi_des(n+1) - phi(n+1)) + K(1,2) * (p_des - p);...
        K(2,1) * (theta_des(n+1) - the(n+1)) + K(2,2) * (q_des - q);...
        K(3,1) * (psi_des(n+1) - psi(n+1)) + K(3,2) * (r_des - r)];
    omega(:,n+1) = omega(:,n) + omega_dot * dt;
    R(:,:,n+1) = R(:,:,n) * expm(Hat(omega(:,n)) * dt);
    v(:,n+1) = v(:,n) + dt*(-m * g * e3 + (T_des * R(:,:,n) * e3)) / m;
    x(:,n+1) = x(:,n) + v(:,n) * dt;
   
end

%% Animation using hgtransform

view(3);
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
set(gca,'XLim',[-5 15],'YLim',[-5 15],'ZLim',[-10 10]);

[x1 y1 z1] = cylinder([0.2 0.2]);
[x2 y2 z2] = cylinder([0.2 0.2]);
g(1) = surface((z1-0.5)*2,y1,x1,'FaceColor','red');
g(2) = surface(x2,(z2-0.5)*2,y2,'FaceColor','blue');

[h1 y3 z3] = cylinder([0.1 0.1]);
[x4 y4 z4] = cylinder([0.1 0.1]);
a(1) = surface(h1-1,y3*2+0.3,z3*0.2+0.1,'FaceColor','red');
a(2) = surface(x4-1,y4*2-0.3,z4*0.2+0.1,'FaceColor','red');
a(3) = surface(h1-1,y3,z3*0.2+0.1,'FaceColor','yellow');

c(1) = surface(h1+1,y3*2+0.3,z3*0.2+0.1,'FaceColor','red');
c(2) = surface(x4+1,y4*2-0.3,z4*0.2+0.1,'FaceColor','red');
c(3) = surface(h1+1,y3,z3*0.2+0.1,'FaceColor','yellow');

b(1) = surface(h1,y3*2+1.3,z3*0.2+0.1,'FaceColor','blue');
b(2) = surface(x4,y4*2+0.7,z4*0.2+0.1,'FaceColor','blue');
b(3) = surface(h1,y3+1,z3*0.2+0.1,'FaceColor','yellow');

d(1) = surface(h1,y3*2-0.7,z3*0.2+0.1,'FaceColor','blue');
d(2) = surface(x4,y4*2-1.3,z4*0.2+0.1,'FaceColor','blue');
d(3) = surface(h1,y3-1,z3*0.2+0.1,'FaceColor','yellow');

quad = hgtransform;
set(g, 'parent', quad);

r1 = hgtransform;
set(a, 'parent', r1);

r2 = hgtransform;
set(b, 'parent', r2);

r3 = hgtransform;
set(c, 'parent', r3);

r4 = hgtransform;
set(d, 'parent', r4);
drawnow

T1x = makehgtform('translate',[-1 0 0]);
T2x = makehgtform('translate',[1 0 0]);

T1y = makehgtform('translate',[0 -1 0]);
T2y = makehgtform('translate',[0 1 0]);

latitude = x(1,:);
longitude = x(2,:);
altitude = x(3,:);
[axis,theta] = axang(R);
ang = [0; 0; 0; 0];

for i = 1:length(latitude)
    translation = makehgtform('translate',latitude(i),longitude(i),altitude(i));
    rotation = makehgtform('axisrotate',axis(:,i),theta(i));
    
    zb = R(:,:,i) * [0;0;1]
    R1 = makehgtform('axisrotate',zb,ang(1));
    R1(1,4) = latitude(i);
    R1(2,4) = longitude(i);
    R1(3,4) = altitude(i);
    
    R2 = makehgtform('axisrotate',zb,ang(2));
    R2(1,4) = latitude(i);
    R2(2,4) = longitude(i);
    R2(3,4) = altitude(i);
    
    R3 = makehgtform('axisrotate',zb,ang(3));
    R3(1,4) = latitude(i);
    R3(2,4) = longitude(i);
    R3(3,4) = altitude(i);
    
    R4 = makehgtform('axisrotate',zb,ang(4));
    R4(1,4) = latitude(i);
    R4(2,4) = longitude(i);
    R4(3,4) = altitude(i);
    
    c1 = R(:,:,i) * [-1; 0; 0];
    c2 = R(:,:,i) * [0; 1; 0];
    c3 = R(:,:,i) * [1; 0; 0];
    c4 = R(:,:,i) * [0; -1; 0];
    
    x1 = makehgtform('translate', [1 0 0]);
    x2 = makehgtform('translate', [-1 0 0]);
    y1 = makehgtform('translate', [0 1 0]);
    y2 = makehgtform('translate', [0 -1 0]);
    
    h1 = makehgtform('translate', [c1(1,1) c1(2,1) c1(3,1)]);
    h2 = makehgtform('translate', [c2(1,1) c2(2,1) c2(3,1)]);
    h3 = makehgtform('translate', [c3(1,1) c3(2,1) c3(3,1)]);
    h4 = makehgtform('translate', [c4(1,1) c4(2,1) c4(3,1)]);
    
    set(quad, 'matrix', translation*rotation);
    set(r1, 'matrix', h1*R1*rotation*x1);
    set(r2, 'matrix', h2*R2*rotation*y2);
    set(r3, 'matrix', h3*R3*rotation*x2);
    set(r4, 'matrix', h4*R4*rotation*y1);
   
    ang = ang + (wp .* dt);
    
    p0 = [latitude(i);longitude(i);altitude(i)];
    p1 = p0 + zb;
    hold on;
    plot3([p0(1), p1(1)], [p0(2), p1(2)], [p0(3), p1(3)]);
    plot3(latitude,longitude,altitude);
    hold on;
    pause(0.01);
end
