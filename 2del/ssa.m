% Inputs:
% c: vector of rate constants; length(c)=m;
% x0: initial conditions; length(x0)=n;
% t_end: end of the time interval [0, t_end]
% t_step: monitoring time step size;
%
% Outputs:
% X: recorded states at these times
%

function [X] = ssa(c,reactants,products,x0,t_end,t_step, e_switch_point, e_switch_d)


[m,n]=size(reactants);

if (m ~= size(products,1) || n~= size(products,2))
    error('Incorrect matrix size of <products>')
end;

if (m ~= length(c))
    error('Incorrect vector size of <c>')
end;

if (n ~= length(x0))
    error('Incorrect vector size of <x0>')
end;

for i=1:n
    if x0(i)<0
        error('Concentrations cannot be negative')
    end;
end;


% initialization:
st = products' - reactants'; % stoichiometry matrix is generated
a = zeros (1,m); % propensity vector
x = x0; % current state
t_proc = 0; % proceeded time
s = floor(t_end / t_step)+1; % no of max. recorded steps
X = zeros (n,s); % records the state
X(:,1) = x0; % state at time zero is x0
i=2; % counts the time steps
m1 = 0; % number of First Order Reactions
m2 = 0; % number of Second Order Reactions
m3 = 0; % number of Homodimer Reactions
m4 = 0; % number of reactions of type 2X + Y -> Z
m5 = 0; % number of reactions of type 2X + 2Y -> Z
rand('state',1);
rand('twister',sum(100*clock))
react_idx = zeros (3,m); % index of one two reactants or threee reactants per reaction
react_type = zeros (1,m); % type of each reaction; 0=FO, 1=SO, 2=HD, 3 = 2X +Y, 4=2X+2Y+Z

% build the two matrices react_id and react_type;
% this also tests for incorrect inputs
for j=1:m % for each reaction
    h1=0; % counts the number of different species involved in reaction j
    h2=0; % counts the total number of reactants involved in reaction j
    
    for k=1:n % for each species
        r = reactants(j,k);
        if r~= 0
            react_idx(h1+1,j)=k;
            h1 = h1+1;
            h2 = h2+r;
        end;        
    end;
       
    if (h1==1 && h2==1)
        react_type(j)=0; %FO
        m1=m1+1;
    elseif (h1==2 && h2==2)
        react_type(j)=1; %SO
        m2=m2+1;
    elseif (h1==1 && h2==2)
        react_type(j)=2; %HD
        m3=m3+1;
	% mmoskon
	elseif (h1==2 && h2==3)
		react_type(j)=3; %2X + Y -> Z
        m4=m4+1;
    elseif (h1==2 && h2==4)
		react_type(j)=4; %2X + 2Y -> Z
        m5=m5+1;     
    elseif (h1==3) && (h2==5)
        react_type(j)=5; %2X + 2Y + Z -> Q
    else
        error('Error');
    end;
end;

c2=[react_type;c]'; % add row and transpose
c2=sortrows(c2,1); % sort rows according to react_type;
c=c2(:,2); % get rid of the first row

% react_idx is a 3xm matrix
react_idx2=[react_type;react_idx]';
react_idx2=sortrows(react_idx2,1);
react_idx=react_idx2(:,(2:4));
react_idx=react_idx';

% st is a nxm matrix
st2=[react_type;st]';
st2=sortrows(st2,1);
st=st2(:,(2:n+1));
st=st';

% auxiliary variables
m1b = m1+1;
m12 = m1+m2;
m12b = m1+m2+1;
m23 = m1+m2+m3;
m23b = m1+m2+m3+1;
m34 = m1+m2+m3+m4;
m34b = m1+m2+m3+m4+1;
m45 = m1+m2+m3+m4+m5;
m45b = m1+m2+m3+m4+m5+1;

% sort the reactants

reactants2 = [react_type;reactants']';
reactants2 = sortrows(reactants2,1);
reactants2 = reactants2(:,(2:end));

a_0 = 0;

for j=1:m1 % for each first order reaction
    a(j) = c(j) * x(react_idx(1,j));
    a_0 = a_0 + a(j);
end;

for j=m1b:m12 % for each second order reaction
    a(j) = c(j) * x(react_idx(1,j)) * x(react_idx(2,j));
    a_0 = a_0 + a(j);
end;

for j=m12b:m23 % for each homodimer
    a(j) = c(j) * x(react_idx(1,j)) * (x(react_idx(1,j))-1) / 2;
    a_0 = a_0 + a(j);
end;

for j=m23b:m34 % for each 2X + Y
    idx1 = react_idx(1,j);
	idx2 = react_idx(2,j);
	
	if (reactants2(j,idx1) == 2)
		reactant1 = x(idx2);
		reactant2 = x(idx1); 	
	else
		reactant1 = x(idx1);
		reactant2 = x(idx2); 	
	end;
	
	a(j) = c(j) * reactant1 * reactant2 * (reactant2 - 1)/2; 
	a_0 = a_0 + a(j);
end;

for j=m34b:m45 % for each 2X + 2Y
    idx1 = react_idx(1,j);
	idx2 = react_idx(2,j);  

    reactant1 = x(idx2);
    reactant2 = x(idx1);

    a(j) = c(j) * reactant1 * (reactant1 - 1)/2 * reactant2 * (reactant2 - 1)/2; 
	a_0 = a_0 + a(j);    
end;


for j=m45b:m % for each 2X + 2Y + Z
    idx1 = react_idx(1,j);
	idx2 = react_idx(2,j);
    idx3 = react_idx(3,j);
    
    if (reactants2(j,idx1) == 1)
		reactant1 = x(idx1);    
		reactant2 = x(idx2); 	
  		reactant3 = x(idx3); 	
    elseif  (reactants2(j,idx2) == 1)
		reactant1 = x(idx2);    
		reactant2 = x(idx1); 	
  		reactant3 = x(idx3); 	
    else
        reactant1 = x(idx3);    
		reactant2 = x(idx1); 	
  		reactant3 = x(idx2); 	
    end;
    
    a(j) = c(j) * reactant1 * reactant2 * (reactant2 - 1) * reactant3 * (reactant3 - 1)/4;
    a_0 = a_0 + a(j);
end;


% this variable saves the state at integer multiples of t_step
x_step = x0;

tt = 0;

while (t_proc < t_end)
    while (tt <= t_step) && (a_0 > 0)
        %determine next step and reaction taking place
        tau = reallog(1/rand)/a_0;
        standard = rand * a_0;
        a_cum = 0; k=0;
        while a_cum < standard
            a_cum = a_cum + a(k+1);
            k=k+1;
        end;
        
        % trigger E injection
        if (t_proc > e_switch_point) && (t_proc < e_switch_point+1)
            x(12) = x(12) + e_switch_d;
        end

        % update state
        x_step = x;
        x = x+st(:,k);

        % update time
        t_proc = t_proc + tau;
        tt = tt + tau;

        % calculate propensities
        a_0 = 0;
        for j=1:m1 % for each first order reaction
            a(j) = c(j) * x(react_idx(1,j));
            a_0 = a_0 + a(j);
        end;

        for j=m1b:m12 % for each second order reaction
            a(j) = c(j) * x(react_idx(1,j)) * x(react_idx(2,j));
            a_0 = a_0 + a(j);
        end;

        for j=m12b:m23 % for each homodimer
            a(j) = c(j) * x(react_idx(1,j)) * (x(react_idx(1,j))-1) / 2;
            a_0 = a_0 + a(j);
        end;
		        
		for j=m23b:m34 % for each 2X + Y
            idx1 = react_idx(1,j);
            idx2 = react_idx(2,j);

            if (reactants2(j,idx1) == 2)
                reactant1 = x(idx2);
                reactant2 = x(idx1); 	
            else
                reactant1 = x(idx1);
                reactant2 = x(idx2); 	
            end;

            a(j) = c(j) * reactant1 * reactant2 * (reactant2 - 1)/2; 
            a_0 = a_0 + a(j);
        end;
        
        for j=m34b:m45 % for each 2X + 2Y
            idx1 = react_idx(1,j);
            idx2 = react_idx(2,j);  

            reactant1 = x(idx2);
            reactant2 = x(idx1);

            a(j) = c(j) * reactant1 * (reactant1 - 1)/2 * reactant2 * (reactant2 - 1)/2; 
            a_0 = a_0 + a(j);    
        end;
        

        for j=m45b:m % for each 2X + 2Y + Z
            idx1 = react_idx(1,j);
            idx2 = react_idx(2,j);
            idx3 = react_idx(3,j);

            if (reactants2(j,idx1) == 1)
                reactant1 = x(idx1);    
                reactant2 = x(idx2); 	
                reactant3 = x(idx3); 	
            elseif  (reactants2(j,idx2) == 1)
                reactant1 = x(idx2);    
                reactant2 = x(idx1); 	
                reactant3 = x(idx3); 	
            else
                reactant1 = x(idx3);    
                reactant2 = x(idx1); 	
                reactant3 = x(idx2); 	
            end;

            a(j) = c(j) * reactant1 * reactant2 * (reactant2 - 1) * reactant3 * (reactant3 - 1)/4;
            a_0 = a_0 + a(j);
        end;
    end;
    
    tt = tt - t_step;
    if a_0==0
        disp('All substrates consumed! Program terminates!')
        disp('Time proceeded so far: ')
        t_proc
        break
    end;

    % record data
    X(:,i) = x_step'; 
 
    i=i+1;
    
end;


for j=i:s
    X(:,j) = x_step;
end;