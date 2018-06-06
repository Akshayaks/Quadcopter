%This function describes the motion of quadrotor as translation motion in
%x y z axis

clc; clear;

Time = 5;
dt = 0.01;
t = 0:dt:Time;
N = length(t);

% constants
m = 2; g = 9.81;
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
R(:,:,1) = vrrotvec2mat([0,1,0,60*pi/180]);
omega(:,1) = zeros(3,1);

T=m*g;
M = zeros(3,N);
e3=[0; 0; 1];

I_inv = inv(I);

for n=1:N-1
    omega(:,n+1) = omega(:,n) + dt*I_inv*(M(:,1)-Hat(omega(:,n))*I*omega(:,n));
    R(:,:,n+1) = R(:,:,n)*expm(Hat(omega(:,n))*dt);
    v(:,n+1) = v(:,n) + dt*(-m*g*e3 + (T*R(:,:,n)*e3))/m;
    x(:,n+1) = x(:,n) + v(:,n)*dt;
end

%% Animation using hgtransform

view(3);
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
set(gca,'XLim',[-5 15],'YLim',[-5 15],'ZLim',[-10 10]);

[x1 y1 z1] = cylinder([0.1 0.1]);
[x2 y2 z2] = cylinder([0.1 0.1]);
h(1) = surface(x1+0.5,y1,z1,'FaceColor','red');
h(2) = surface(z2,x2,y2+0.5,'FaceColor','blue');

combinedobject = hgtransform;
set(h, 'parent', combinedobject)
drawnow

latitude = x(1,:);
longitude = x(2,:);
altitude = x(3,:);
[axis,theta] = axang(R);

for i = 1:length(latitude)
    
    translation = makehgtform('translate',[latitude(i),longitude(i),altitude(i)]);
    rotation = makehgtform('axisrotate',axis(:,i),theta(i));
    set(combinedobject, 'matrix', translation*rotation);
    pause(0.1);
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
