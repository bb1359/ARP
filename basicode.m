clc
clear

ode_start = 0;
ode_end = 60;
ode_n = 5000;
e_switch_point = 20;
e_switch_d = 1000;

t2 = linspace(ode_start,ode_end,ode_n)';
frst = 1;

if(frst)
  startValues = [1 0.3 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
  fun = @(x, t) moskon_full(x,t,e_switch_point,e_switch_d);
  toplot = [13 14 15 16];
else
  startValues = zeros(10,1);
  fun = @(x, t) moskon_simple(x,t,e_switch_point,e_switch_d);
  toplot = [4 5 6 7];
end

tic
if(exist('OCTAVE_VERSION', 'builtin') ~= 0)
  x2 = lsode (fun, startValues, t2);
  graphics_toolkit("gnuplot")
else
  x2 = ode45(fun, [ode_start;ode_end;ode_n], startValues);
end
toc

p = plot(t2,x2(:, toplot));
legend('IFN', 'A2', 'B2', 'E2');


