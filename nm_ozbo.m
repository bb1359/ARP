function [iszero, solution, iter] = nm_ozbo(initial_solution, fobj, fden)
  % initialise n+1 solutions
  initials = ((rand(length(initial_solution))-0.5)*0.01) + initial_solution;
  initials = [initial_solution; initials];

  % define function for nelder_mead
  function [sol,iter] = fobjden(sol) 
    sol = fden(sol);
    [sol,iter] = fobj(sol);
    #disp(err);
    #fflush(stdout);
  end;

  % start nelder_mead algorithm
  f = @(x)fobjden(x);
  [solution, iter] = nelder_mead(initials, f);
  
  if (iter > 0)
    iszero = 1;
  else
    iszero = 0;
  end;

end;

  
