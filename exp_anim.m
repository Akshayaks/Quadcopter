%This function describes the motion of quadrotor as translation motion in
%x y z axis

clc; clear;

Time = 5;
dt = 0.01;
t = 0:dt:Time;
N = length(t);

% constants
m = 2; 
g = 9.81;
L = 0.5;
I = diag([0.081, 0.0812, 0.1320]);
% b = 3.2*10^-6;
Jr = 5 * 10^-5;

%of the blade
area = 10^-3; 
radius = 0.03;
L = 0.08; %distance from the axis of rotation to the rotors

%Aerodynamic constants
rho = 1;
Cd = 0.19; Ct = 0.23;


% F = Kf*omega*omega;
% M = Km*omega*omega;
% M = matrix * omega;

Kf = 0.5 * rho * area * Ct * radius * radius;
Km = 0.5 * rho * area * Cd * radius * radius;
matrix = [Kf Kf Kf Kf;0 Kf*L 0 -Kf*L;-Kf*L 0 Kf*L 0;Km -Km Km -Km];

% state
x = zeros(3,N);
v = zeros(3,N);
R = zeros(3,3,N);
omega = zeros(3,N);
w = [3100;3100;3100;4718.62];

% initialize
x(:,1) =[5; 5; 5];
v(:,1)=[5; 0; 0];
%R(:,:,1) = [1 0 0;0 1 0;0 0 1];
R(:,:,1) = vrrotvec2mat([0,1,0,60*pi/180]);
omega(:,1) = zeros(3,1);

T=m*g;
M = zeros(3,N);
e3=[0; 0; 1];

I_inv = inv(I);

for n=1:N-1
    V = matrix * w .^ 2;
    T = V(1);
    M = V(2:end);
    gr = Jr * (Hat(omega(:,n)) * e3) * (w(1) - w(2) + w(3) - w(4));
    omega(:,n+1) = omega(:,n) + dt * I_inv * (M(:,1) - Hat(omega(:,n)) * I * omega(:,n) + gr);
    R(:,:,n+1) = R(:,:,n) * expm(Hat(omega(:,n)) * dt);
    v(:,n+1) = v(:,n) + dt*(-m * g * e3 + (T * R(:,:,n) * e3)) / m;
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
ang = 20;
wr = 20;

for i = 1:length(latitude)
    translation = makehgtform('translate',latitude(i),longitude(i),altitude(i));
    rotation = makehgtform('axisrotate',axis(:,i),theta(i));
    
    zb = R(:,:,i) * [0;0;1];
    R1 = makehgtform('axisrotate',zb,ang);
    set(r1, 'matrix', R1);
    R1(1,4) = latitude(i);
    R1(2,4) = longitude(i);
    R1(3,4) = altitude(i);
    
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
    set(r2, 'matrix', h2*R1*rotation*y2);
    set(r3, 'matrix', h3*R1*rotation*x2);
    set(r4, 'matrix', h4*R1*rotation*y1);
    
    ang = ang + wr * 0.3;
    
    p0 = [latitude(i);longitude(i);altitude(i)];
    p1 = p0 + axis(:,i);
    hold on;
    plot3([p0(1), p1(1)], [p0(2), p1(2)], [p0(3), p1(3)]);
    plot3(latitude,longitude,altitude);
    hold on;
    pause(0.3);
end

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

% plot(t,x);
% grid on;
% title('x')
% xlabel('time (s)'); ylabel('x (m)');
% 
% figure
% plot(t,v);
% grid on;
% title('v')
% xlabel('time (s)'); ylabel('v (m/s)');
% 
% figure
% plot(t,omega);
% grid on;
% title('\omega')
% xlabel('time (s)'); ylabel('\omega (rad/s)');
