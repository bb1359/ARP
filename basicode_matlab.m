%clc
%clear

% x =
% [1  2  3  4    5    6    7    8    9   ...]
% [pA pB pI pAA2 pAB2 pBA2 pBE2 pIB2 pII ...]
% [... 10 11 12 13 14 15 16 17 18 19]
% [... A  B  E  I  A2 B2 E2 mA mB mI]

function dx = moskon(x,t)
    pA = x(1); pB = x(2); pI = x(3);
    pAA2 = x(4); pAB2 = x(5);
    pBA2 = x(6); pBE2 = x(7);
    pIB2 = x(8); pII = x(9);
    A = x(10); B = x(11); E = x(12); I = x(13);
    A2 = x(14); B2 = x(15); E2 = x(16);
    mA = x(17); mB = x(18); mI = x(19);
    
    k1=1;k2=1;k3=1;k4=1;k5=1;k6=1;
    k_1=1;k_2=1;k_3=1;k_5=1;k_6=1;
    k_4=-1;
    k7=1e-3;k8=k7;k9=k7;
    k_7=1;k_8=1;k_9=1;
    b_a=10;b_b=30;b_i=10;
    fi=50;fi_a=50;
    a_a=5;a_b=1;a_i=5;
    g_a=10;g_b=10;g_i=10;
    d_a=1;d_b=1;d_i=1;d_e=1;
    
    dx(1) = -k1*pA*A2 + k_1*pBA2 - k2*pA*B2 + k_2*pAB2;
    dx(2) = -k3*pB*A2 + k_3*pBA2 - k4*pB*E2 + k_4*pBE2;
    dx(3) = -k5*pI*B2 + k_5*pIB2 - k6*pI*I  + k_6*pII;
    
    dx(4) = k1*pA*A2 - k_1*pAA2;
    dx(5) = k2*pA*B2 - k_2*pAB2;
    dx(6) = k3*pB*A2 - k_3*pBA2;
    dx(7) = k4*pB*E2 - k_4*pBE2;
    dx(8) = k5*pI*B2 - k_5*pIB2;
    dx(9) = k6*pI*I  - k_6*pII;
    
    dx(10) = -2*k7*A^2 + 2*k_7*A2 + g_a*mA - d_a*A;
    dx(11) = -2*k8*B^2 + 2*k_8*B2 + g_b*mB - d_b*B;
    dx(12) = -2*k9*E^2 + 2*k_9*E2          - d_e*E;
    dx(13) = -k6*I*pI  + k_6*pII  + g_i*mI - d_i*I;
    
    dx(14) = -k1*pA*A2 + k_1*pAA2 - k3*B2*pA  + k_3*pBA2 + k7*A^2 - k_7*A2 - d_a*A2;
    dx(15) = -k2*pA*B2 + k_2*pAB2 - k5*B2*pI  + k_5*pIB2 + k8*B^2 - k_8*B2 - d_b*B2;
    dx(16) = -k4*pB*E2 + k_4*pBE2                        + k9*E^2 - k_9*E2 - d_e*E2;
    
    dx(17) = b_a*pA + b_a*fi_a*pAA2 - a_a*mA;
    dx(18) = b_b*pB + b_b*fi  *pBE2 - a_b*mB;
    dx(19) = b_i*pI + b_i*fi  *pIB2 - a_i*mI;
end

function xdot = f (x,t)
    r = 0.25;
    k = 1.4;
    a = 1.5;
    b = 0.16;
    c = 0.9;
    d = 0.8;
    xdot(1) = r*x(1)*(1 - x(1)/k) - a*x(1)*x(2)/(1 + b*x(1));
    xdot(2) = c*a*x(1)*x(2)/(1 + b*x(1)) - d*x(2);     
end
 
%
% Using LSODE again, with initial conditions x(1)=1 and x(2)=2 on [0,50] with 200 points:
%

%OCTAVE implementation
%x = lsode ("f", [1; 2], (t = linspace (0, 50, 200)'));
%
%MATLAB implementation
function izracunaj
    ode45(@f, [1;2], [0;50;200]);
    %ode45(@moskon, ones(19,1), [0;1;1000]');
end

%x2 = lsode ("moskon", ones(19,1), (t2=linspace(0,1,1000)'));

%graphics_toolkit("gnuplot")
%#set term dumb;
%plot(t2,x2);

