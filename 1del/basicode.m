clc
clear

ode_start = 0;
ode_end = 100;
ode_n = 500;

e_switch_point = 20;
e_switch_d = 1000;
pA_init=1;pB_init=0.4;pI_init=1;
k1=1;k2=1;k3=1;k4=1;k5=1;k6=1;
k_1=1;k_2=1;k_3=1;k_5=1;k_6=1;
k_4=0.1;
k7=1e-3;k8=k7;k9=k7;
k_7=1;k_8=1;k_9=1;
b_a=10;b_b=30;b_i=10;
fi=50;fi_a=50;
a_a=5;a_b=1;a_i=5;
g_a=10;g_b=10;g_i=10;
d_a=1;d_b=1;d_i=1;d_e=1;

args = [e_switch_point;e_switch_d;pA_init;pB_init;pI_init;k1;k2;k3;k4;k5;k6;
        k_1;k_2;k_3;k_5;k_6;k_4;k7;k8;k9;k_7;k_8;k_9;b_a;b_b;b_i;fi;fi_a;
        a_a;a_b;a_i;g_a;g_b;g_i;d_a;d_b;d_i;d_e]';

t2 = linspace(ode_start,ode_end,ode_n)';

% run full
%what2run = 0;
% run simple 
what2run = 1;
% run simpler
%what2run = 2;
% run simplest
%what2run = 3;

if(what2run == 0)
  startValues = [pA_init pB_init pI_init 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
  fun = @(x, t) moskon_full(x,t,args);
  toplot = [13 14 15 16];
elseif(what2run == 1)
  startValues = zeros(10,1);
  fun = @(x, t) moskon_simple(x,t,args);
  toplot = [4 5 6 7];
else
  exit(1);
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
