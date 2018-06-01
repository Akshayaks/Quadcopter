%This function is for only translation motion along z-axis 
%using ode45 function

function hover(mass)
  g = 10;
  %rot_mat = calc
  rot_mat = [1 0 0;0 1 0;0 0 1];
  v = rot_mat * [0;0;1];
  v0 = 0;
  T = mass * g;
  odefunv = @(t) g - (T * v(3,1))/mass;
  t_span = [0:100];
  [tSol vSol] = ode45(odefunv,t_span,v0);
  plot(tSol,vSol);
  hold on;

   s = size(vSol);
   for i = 1:s %As v changes x needs to be calculated every time
     odefunx = @(t) vSol(i);
     x0 = 10;
     [t xSol] = ode45(odefunx,t_span,x0);
     axis([0 10 0 20]);
     x_value(i) = xSol(i);
   endfor
   plot(tSol,x_value);
endfunction

function R = calc
  A = rand(3);
  R = A - A';
endfunction
