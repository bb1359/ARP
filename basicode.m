clc
clear

ode_start = 0;
ode_end = 100;
ode_n = 1000;

t2 = linspace(ode_start,ode_end,ode_n)';
startValues = [1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';

if(exist('OCTAVE_VERSION', 'builtin') ~= 0)
  lsode_options('maximum step size', 0.99);
  x2 = lsode ("moskon_full", startValues, t2);
  graphics_toolkit("gnuplot")
else
  x2 = ode45(@moskon_full, [ode_start;ode_end;ode_n], startValues);
end

p = plot(t2,x2(:, [14 15 16 13]));
legend('A2', 'B2', 'E2', 'IFN');
% legend(p, "pA", "pB", "pI", 
%           "pAA2", "pAB2", "pBA2", "pBE2", "pIB2", "pII", 
%           "A",  "B",  "E",  "I",  
%           "A2", "B2", "E2", 
%           "mA", "mB", "mI");
