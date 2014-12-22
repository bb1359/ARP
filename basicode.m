clc
clear

ode_start = 0;
ode_end = 50;
ode_n = 5000;

t2 = linspace(ode_start,ode_end,ode_n)';
frst = 0;

if(frst)
  startValues = [1 0.3 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
  fun = "moskon_full";
  fun_m = @moskon_full;
  toplot = [13 14 15 16];
else
  startValues = zeros(10,1);
  fun = "moskon_simple";
  fun_m = @moskon_simple;
  toplot = [4 5 6 7];
end


if(exist('OCTAVE_VERSION', 'builtin') ~= 0)
  x2 = lsode (fun, startValues, t2);
  graphics_toolkit("gnuplot")
else
  x2 = ode45(fun_m, [ode_start;ode_end;ode_n], startValues);
end

p = plot(t2,x2(:, toplot));
legend('IFN', 'A2', 'B2', 'E2');


