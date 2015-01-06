% stochastic model for ebola treatment
% [1  2  3  4    5    6    7    8    9   ...]
% [pA pB pI pAA2 pAB2 pBA2 pBE2 pIB2 pII ...]
% [... 10 11 12 13 14 15 16 17 18 19]
% [... A  B  E  I  A2 B2 E2 mA mB mI]
% Inputs:
% t_end: end of the time interval [0, t_end] (period of simulation)
% t_step: monitoring time step size;
% X_init: initial conditions
% observ: observed chemical species
% plot_on: turns the plotting on
%
% Outputs:
% X_mean: time course of observed chemical species
%
% sample call
% inverter_ssa(10000,10,[50 1 0 0 0]',[5],1);
% inverter_ssa(10000,10,[600 1 0 25 1200]',[5],1);
%function [X_mean] = inverter_ssa(t_end,t_step,X_init,observ,plot_on)

pA_init=1;pB_init=0.4;pI_init=1;
t_end = 10000;
t_step = 10;
% IFN, A2, B2, E2
observ = 13:16;
plot_on = true;

% average of how many runs?
runs = 1;

disp('SSA - ebola');

% reaction volume (omega = AV)
A=6.0221415*10^23;
V=2*10^(-15); % cell volume
OMEGA=A*V*10^(-9);  

% parameters
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

% -> define the reactions (i.e. the reactants and the products matrix)
% [1  2  3  4    5    6    7    8    9   ...]
% [pA pB pI pAA2 pAB2 pBA2 pBE2 pIB2 pII ...]
% [... 10 11 12 13 14 15 16 17 18 19]
% [... A  B  E  I  A2 B2 E2 mA mB mI]
% base for vector construction
base_vector = zeros(1,19);
% rx are reactants, px are products or reactions 
% generate c from k and type of reaction

% pA + A2 => pAA2
r1 = base_vector;   p1 = base_vector;
r1([1,14]) = 1;     p1(4) = 1;
c1 = k1/OMEGA;
% pAA2 => pA + A2
r2 = p1;            p2 = r1;
c2 = k_1;

% pA + B2 => pAB2
r3 = base_vector;   p3 = base_vector;
r3([1,15]) = 1;     p3(5) = 1;
c3 = k2/OMEGA;
% pAB2 => pA + B2
r4 = p3;            p4 = r3;
c4 = k_2;

% pB + A2 => pBA2
r5 = base_vector;   p5 = base_vector;
r5([2,14]) = 1;     p5(6) = 1;
c5 = k3/OMEGA;
% pBA2 => pB + A2
r6 = p5;            p6 = r5;
c6 = k_3;

% pB + E2 => pBE2
r7 = base_vector;   p7 = base_vector;
r7([2,16]) = 1;     p7(7) = 1;
c7 = k4/OMEGA;
% pBE2 => pB + E2
r8 = p7;            p8 = r7;
c8 = k_4;

% pI + B2 => pIB2
r9 = base_vector;   p9 = base_vector;
r9([3,15]) = 1;     p9(8) = 1;
c9 = k5/OMEGA;
% pIB2 => pI + B2
r10 = p9;           p10 = r9;
c10 = k_5;

% pI + I => pII
r11 = base_vector;  p11 = base_vector;
r11([3,13]) = 1;    p11(9) = 1;
c11 = k6/OMEGA;
% pII => pI + I
r12 = p11;          p12 = r11;
c12 = k_6;

% 2A => A2
r13 = base_vector;  p13 = base_vector;
r13(10) = 2;        p13(14) = 1;
c13 = 2*k7/OMEGA;
% A2 => 2A
r14 = p13;          p14 = r13;
c14 = k_7;

% 2B => B2
r15 = base_vector;  p15 = base_vector;
r15(11) = 2;        p15(15) = 1;
c15 = 2*k8/OMEGA;
% B2 => 2B
r16 = p15;          p16 = r15;
c16 = k_8;

% 2E => E2
r17 = base_vector;  p17 = base_vector;
r17(12) = 2;        p17(16) = 1;
c17 = 2*k9/OMEGA;
% E2 => 2E
r18 = p17;          p18 = r17;
c18 = k_9;

% pA => pA + mA
r19 = base_vector;  p19 = base_vector;
r19(1) = 1;         p19([1,17]) = 1;
c19 = b_a;

% pAA2 => pAA2 + mA
r20 = base_vector;  p20 = base_vector;
r20(4) = 1;         p20([4,17]) = 1;
c20 = b_a * fi_a;

% pB => pB + mB
r21 = base_vector;  p21 = base_vector;
r21(2) = 1;         p21([2,18]) = 1;
c21 = b_b;

% pBE2 => pBE2 + mB
r22 = base_vector;  p22 = base_vector;
r22(7) = 1;         p22([7,18]) = 1;
c22 = b_b * fi;

% pI => pI + mI
r23 = base_vector;  p23 = base_vector;
r23(3) = 1;         p23([3,19]) = 1;
c23 = b_i;

% pIB2 => pIB2 + mI
r24 = base_vector;  p24 = base_vector;
r24(8) = 1;         p24([8,19]) = 1;
c24 = b_i * fi;

% mA => mA + A
r25 = base_vector;  p25 = base_vector;
r25(17) = 1;        p25([10,17]) = 1;
c25 = g_a;

% mB => mB + B
r26 = base_vector;  p26 = base_vector;
r26(18) = 1;        p26([11,18]) = 1;
c26 = g_b;

% mI => mI + I
r27 = base_vector;  p27 = base_vector;
r27(19) = 1;        p27([13,19]) = 1;
c27 = g_i;

% A => 0
r28 = base_vector;  p28 = base_vector;
r28(10) = 1;
c28 = d_a;

% A2 => 0
r29 = base_vector;  p29 = base_vector;
r29(14) = 1;
c29 = d_a;

% B => 0
r30 = base_vector;  p30 = base_vector;
r30(11) = 1;
c30 = d_b;

% B2 => 0
r31 = base_vector;  p31 = base_vector;
r31(15) = 1;
c31 = d_b;

% E => 0
r32 = base_vector;  p32 = base_vector;
r32(12) = 1;
c32 = d_e;

% E2 => 0
r33 = base_vector;  p33 = base_vector;
r33(16) = 1;
c33 = d_e;

% I => 0
r34 = base_vector;  p34 = base_vector;
r34(13) = 1;
c34 = d_i;

% mA => 0
r35 = base_vector;  p35 = base_vector;
r35(17) = 1;
c35 = a_a;

% mB => 0
r36 = base_vector;  p36 = base_vector;
r36(18) = 1;
c36 = a_b;

% mI => 0
r37 = base_vector;  p37 = base_vector;
r37(19) = 1;
c37 = a_i;

reactants = [r1;r2;r3;r4;r5;
            r6;r7;r8;r9;r10;
            r11;r12;r13;r14;r15;
            r16;r17;r18;r19;r20;
            r21;r22;r23;r24;r25;
            r26;r27;r28;r29;r30;
            r31;r32;r33;r34;r35;
            r36;r37];

products = [p1;p2;p3;p4;p5;
            p6;p7;p8;p9;p10;
            p11;p12;p13;p14;p15;
            p16;p17;p18;p19;p20;
            p21;p22;p23;p24;p25;
            p26;p27;p28;p29;p30;
            p31;p32;p33;p34;p35;
            p36;p37];

constants = [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20, c21, c22, c23, c24, c25, c26, c27, c28, c29, c30,c31, c32, c33, c34, c35, c36, c37];

% initialization
X_init = base_vector;
X_init(1) = pA_init;
X_init(2) = pB_init;
X_init(3) = pI_init;

s = floor(t_end / t_step)+1; % no of max. recorded steps
t = 0:t_step:(s-1)*t_step; % sample points are predefined by the step_size
X_total = zeros(runs,length(observ),s); % data storage for all observed species over all runs

for k=1:runs
    [X] = ssa(constants,reactants,products,X_init,t_end,t_step);
   
    for l=1:length(observ)
        X_total(k,l,:)=X(observ(l),:);
    end;
end;

X_mean = mean(X_total,1);
X_means = squeeze(X_mean);

if (plot_on==1)

    
    figure(1), plot(t, X_means(:,:));
    xlabel('Time [s]');
    ylabel('Number of molecules');
    box off;
    
end