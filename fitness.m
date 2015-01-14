% function that returns fitness function with given 
% starting parameters

function [err y] = fitness(initargs,fun,start_values,args)
  [err y] = fitness_help(initargs,fun,start_values,args);
  if (err < 0)
    err = 999999999;
  end
end

function [err y] = fitness_help(initargs,fun,start_values,args)
  y = 0;

  initargsCell = num2cell(initargs);
  [ode_start ode_end ode_n i_min i_max e_time ifn saving]=initargsCell{:};

  tt = linspace(ode_start,ode_end,ode_n);
  fun = @(x,t) fun (x,t,args);
  [xx istat msg] = lsode (fun, start_values, tt);

  e7 = 10000;
  %error at solving
  if (istat != 2)
    err = -1;
    return;
  end
  y = xx(:, [ifn:(ifn+3)]);
  xx = xx(:, ifn);

  %Did it actually reach i_min?
  reach_i_min  = find(diff( xx>i_min ));
  if(size(reach_i_min)(1) == 0)
    err = 10*e7 + sum(i_min - xx);
    return;
  end;

  reach_i_min = reach_i_min(1);

  %%Did it went 4* over i_max?
  %if (max(xx) > (i_max*4))
  %  err = 8*e7 + sum(abs(xx - i_max));
  %  return;
  %else
  %end

  % raise if IFN already high at switch
  e_idx = find(diff(tt > e_time));
  if (xx(e_idx) > i_min/2) 
    err = 5*e7 + sum(abs(i_max-xx));
    return;
  end

  % if ifn raises before e_signal or
  % it raise 15 hours or more after...
  if ((tt(reach_i_min) > e_time+15) || (tt(reach_i_min) < e_time))
    err = 5*e7 + sum(abs(i_max - xx));
    return;
  end

  xx = xx(reach_i_min(1)+1:end);
  err = [xx(find(xx>i_max))-i_max; xx(find(xx<i_min))-i_min];
  err = (1/ode_n) * sqrt(sum(err.^2));
end

