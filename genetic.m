function [pars,best,fmin] = genetic(fobj,fden,flocal,nd,Tol)
  if nargin != 5,
    disp('Number of args: 5');
    return;
  end

  fmin = Tol+1;
  iter = 0;
  II = 0;
  new_candidates = generate_perms(nd);
  pars = [];

  while(fmin > Tol)
    fall = [];

    for c = new_candidates';
      denorm_args = fden(c');
      [fnew y] = fobj(denorm_args);
      fall = [fall; fnew];
      %printf('\t\t%d: %10.5f\n', iter, fnew);
      %fflush(stdout);

      if(fnew < Tol)
        pars = c';
        fmin = fnew;
        best = y;
        break;
      end;
      iter = iter + 1;
    end
    if(fnew < Tol)
      break;
    end;

    %find best 5 amongst new_candidates
    best = find(sort(fall)(6) > fall);

    %pair and mutate best
    if(iter < 200)
      new_candidates = pair_and_mutate(new_candidates(best,:), 33, 0.01);
    else
      II = II + iter;
      iter = 0;
      new_candidates = generate_perms(nd);
    end;
  end

  iter = iter + II;
  printf('iter: %3d\n', iter);
end;

function perms = generate_perms(n)
  perms = [];
  for _ = 1:n
    perms = [perms; randperm(n)-1];
  end
  perms = (perms' + rand(n))/n;
end;

function new_cand = pair_and_mutate(old, n, mut)
  new_cand = [];
  sz = size(old);
  for _ = 1:n
    choose = [floor(rand(1,sz(2)) * sz(1)) + 1; 1:sz(2)];
    mutate =  (rand(1,sz(2)) - 0.5) * mut;

    %create mask
    mask = zeros(sz);
    for c = choose
      mask( c(1), c(2)) = 1;
    end

    %apply mask and mutate
    new_cand = [new_cand; sum(old.*mask) + mutate];
  end

  new_cand = simplebounds(new_cand);
end;

function new_s = simplebounds(s)
  % Apply the lower bound
  new_s=s;
  I=new_s<0;
  new_s(I)=0;
  
  % Apply the upper bounds 
  J=new_s>1;
  new_s(J)=1;
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

  
