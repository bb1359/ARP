function [bestnest,fmin]=cuckoo_search(n,fobj,fden,nd,pa,Tol)
  if nargin != 6,
    disp('Number of args: 6');
    return;
  end

  % Random initial solutions
  for i=1:n,
    nest(i,:)=rand(nd,1);
  end

  % Get the current best
  fitness=10^10*ones(n,1);
  [fmin,bestnest,nest,fitness]=get_best_nest(nest,nest,fitness,fobj,fden);

  N_iter=n;
  %% Starting iterations
  while (fmin>Tol),
  
    % Generate new solutions (but keep the current best)
    new_nest=get_cuckoos(nest,bestnest);
    [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness,fobj,fden);
    % Update the counter
    N_iter=N_iter+n; 
    % Discovery and randomization
    new_nest=empty_nests(nest,pa) ;
    
    % Evaluate this set of solutions
    [fnew,best,nest,fitness]=get_best_nest(nest,new_nest,fitness,fobj,fden);
    % Update the counter again
    N_iter=N_iter+n;
    % Find the best objective so far  
    if fnew<fmin,
      fmin=fnew;
      bestnest=best;
    end

    %printf('... Iterations %2d\n', N_iter/(2*n));
    %fflush(stdout);
    
  end %% End of iterations

  %% Post-optimization processing
  %% Display all the nests
  
  name = sprintf('%016.5f', fmin);
  disp(name);
  printf('iterations: %3d\n', N_iter);
  fflush(stdout);
  %disp(strcat('iterations: ',num2str(N_iter)));

%% Get cuckoos by ramdom walk
function nest=get_cuckoos(nest,best)
  % Levy flights
  n=size(nest,1);
  % Levy exponent and coefficient
  % For details, see equation (2.21), Page 16 (chapter 2) of the book
  % X. S. Yang, Nature-Inspired Metaheuristic Algorithms, 2nd Edition, Luniver Press, (2010).
  beta=3/2;
  sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
  
  for j=1:n,
      s=nest(j,:);
      
      % This is a simple way of implementing Levy flights
      % For standard random walks, use step=1;
      %% Levy flights by Mantegna's algorithm
      u=randn(size(s))*sigma;
      v=randn(size(s));
      step=u./abs(v).^(1/beta);
    
      % In the next equation, the difference factor (s-best) means that 
      % when the solution is the best solution, it remains unchanged.     
      stepsize=0.01*step.*(s-best);
      % Here the factor 0.01 comes from the fact that L/100 should the typical
      % step size of walks/flights where L is the typical lenghtscale; 
      % otherwise, Levy flights may become too aggresive/efficient, 
      % which makes new solutions (even) jump out side of the design domain 
      % (and thus wasting evaluations).
      
      % Now the actual random walks or flights
      s=s+stepsize.*randn(size(s));

      % Apply simple bounds/limits
      nest(j,:)=simplebounds(s);
  end

%% Find the current best nest
function [fmin,best,nest,fitness]=get_best_nest(nest,newnest,fitness,fobj,fden)
  % Evaluating all new solutions
  for j=1:size(nest,1),
      denorm_args=fden(newnest(j,:));
      fflush(stdout);
      [fnew y] = fobj(denorm_args);

      %name = sprintf('%011.5f',fnew);
      %disp(name);

      %if(fnew < 1000)
      %  plot(y);
      %  axis([0 size(y,1) 0 2000]);
      %  legend('IFN','A2','B2','E2');
      %  %set(get(gca,'YLabel'),'Rotation',0)
      %  %ylabel(num2cell(denorm_args));
      %  args = newnest(j,:);
      %  xlabel({num2str(args(1:11)), num2str(args(12:22)), num2str(args(23:end))}, 'FontSize',5);

      %  print (strcat('md5imgs/', name, '_', (num2hex(rand))(13:end), '.svg'));
      %end

      if fnew<=fitness(j),
         fitness(j)=fnew;
         nest(j,:)=newnest(j,:);
      end
  end
  % Find the current best
  [fmin,K]=min(fitness) ;
  best=nest(K,:);

%% Replace some nests by constructing new solutions/nests
function new_nest=empty_nests(nest,pa)
  % A fraction of worse nests are discovered with a probability pa
  n=size(nest,1);
  % Discovered or not -- a status vector
  K=rand(size(nest))>pa;

  % In the real world, if a cuckoo's egg is very similar to a host's eggs, then 
  % this cuckoo's egg is less likely to be discovered, thus the fitness should 
  % be related to the difference in solutions.  Therefore, it is a good idea 
  % to do a random walk in a biased way with some random step sizes.  
  %% New solution by biased/selective random walks
  stepsize=rand*(nest(randperm(n),:)-nest(randperm(n),:));
  new_nest=nest+stepsize.*K;
  for j=1:size(new_nest,1)
      s=new_nest(j,:);
    new_nest(j,:)=simplebounds(s);  
  end

% Application of simple constraints
function s=simplebounds(s)
  % Apply the lower bound
  ns_tmp=s;
  I=ns_tmp<0;
  ns_tmp(I)=0;
  
  % Apply the upper bounds 
  J=ns_tmp>1;
  ns_tmp(J)=1;
  % Update this new move 
  s=ns_tmp;

