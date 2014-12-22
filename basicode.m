clc
clear

ode_start = 0;
ode_end = 1000;
ode_n = 200000;

t2=linspace(ode_start,ode_end,ode_n)';

if(exist('OCTAVE_VERSION', 'builtin') ~= 0)
  x2 = lsode ("moskon_full", ones(19,1), t2);
  graphics_toolkit("gnuplot")
else
  x2 = ode45(@moskon_full, [ode_start;ode_end;ode_n], ones(19,1));
end

p = plot(t2,x2);
legend(p, "pA", "pB", "pI", 
          "pAA2", "pAB2", "pBA2", "pBE2", "pIB2", "pII", 
          "A",  "B",  "E",  "I",  
          "A2", "B2", "E2", 
          "mA", "mB", "mI");

