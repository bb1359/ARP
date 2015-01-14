function [best,fmin] = hybercube(fobj,fden,flocal,nd,Tol)
  if nargin != 5,
    disp('Number of args: 5');
    return;
  end

  fmin = Tol+1;
  iter = 0;

  while(fmin > Tol)
    new_candidates = generate_perms(nd);
    for c = new_candidates;
      denorm_args = fden(c');
      [fnew y] = fobj(denorm_args);

      if(fnew < Tol)
        fmin = fnew;
        name = sprintf('%016.5f',fnew);
        disp(name);
        printf('iterations: %3d\n', iter);
        fflush(stdout);
        break;
      end;
      iter = iter + 1;
    end
    
    if(fmin < Tol)
      break;
    end;
  end

  best = y;
  return;
end;

function perms = generate_perms(n)
  perms = [];
  for _ = 1:n
    perms = [perms; randperm(n)];
  end
  perms = perms';

  perms = perms + rand(n);
  perms = perms/n;
end;

function perms = generate_sudoku_perms(cols,rows)
  perms = zeros(rows,cols);

  for r = 1:rows
    flag = 1;
    while(flag)
      flag = 0;
      R = randperm(cols);
      for testRow = 1:r
        flag = flag + sum(perms(testRow,:) == R);
      end;
      if (~flag)
        perms(r,:) = R;
      end;
    end;
  end;

  perms = perms/cols;
end
