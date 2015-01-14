function [iszero, solution, iter] = seminar3(combo)

format long
graphics_toolkit('gnuplot')
lsode_options('minimum step size', 1e-13);
lsode_options('maximum step size', 0.05);
lsode_options('step limit', 10000);

save_history = true;

ode_start = 0;
ode_end = 100;
ode_n = 500;

i_min = 250;
i_max = 500;

e_start = 20;
e_len = 1;
e_d = 500;

pA_init=1;
pB_init=1;
pI_init=1;

initargs = [e_start e_len e_d pA_init pB_init pI_init];

what2run = 2;

%%%%%%%%%%%%%%%%%%

if(what2run == 0)
  start_values = [pA_init pB_init pI_init zeros(1,16)]';
  fun = @(x,t,args)moskon_full(x,t,args,initargs);
  ifn = 13; %[13 14 15 16];
elseif(what2run == 1)
  start_values = zeros(10,1);
  fun = @(x,t,args)moskon_simple(x,t,args,initargs);
  ifn = 4; %[4 5 6 7];
elseif(what2run == 2)
  start_values = zeros(7,1);
  fun = @(x,t,args)moskon_simpler(x,t,args,initargs);
  ifn = 1; 
else
  exit(1);
end

denorm_exp = [8;8;8;8;8;8;8;8;8;8;8;8;8;8;8;8;8;8;
              9;9;9;3;3;6;6;6;9;9;9;6;6;6;6]';
denorm_off = [-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;-2;
              -6;-6;-6;0;0;-3;-3;-3;-6;-6;-6;-3;-3;-3;-3]';

fitness_func = @(args) fitness([ode_start ode_end ode_n i_min i_max, e_start, ifn, save_history], fun, start_values, args);
denorm_func = @(x) denorm(x,denorm_exp,denorm_off);

%args = [0.25000;0.25000;0.25000;0.25000;0.25000;0.25000;0.25000;0.25000;0.25000;0.25000;0.25000;0.12500 ;-0.12500;-0.12500;-0.12500;0.25000;0.25000;0.25000;0.77778;0.83079;0.77778;0.56632;0.56632;0.61650;0.50000;0.61650;0.77778;0.77778;0.77778;0.50000;0.50000;0.50000;0.50000]';

% Number of nests
n = 25;

% Discovery rate of alien eggs/solutions
pa=0.25;

%% Change this if you want to get better results
% Tolerance
Tol=1000;

% cuckoo search
%for _ = 1:1000
%  [bestnest,fmin] = cuckoo_search(n,fitness_func,denorm_func,33,pa,Tol);
%end;
%[err y] = fitness_func(denorm_func(bestnest));

%hybnercube global search
%for b = 1:1000
%  [best, fmin] = hybercube(fitness_func, denorm_func, 1, 33, 1000);
%end;

%genetics
%for b = 1:1000
  %[pars, best, fmin] = genetic(fitness_func, denorm_func, 1, 33, 100);
  %printf('pars '); disp(pars);
  %printf('\nerr %2.4f\n',fmin);
  %fflush(stdout);
%end;

%nelder mead
[iszero, solution, iter] = nm_ozbo(combo, fitness_func, denorm_func);


%[err y] = fitness_func(denorm_func(solution));
%plot(y);
%legend('IFN','A2','B2','E2');
end;

