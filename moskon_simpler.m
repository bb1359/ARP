% x =
% [1  2  3  4  5  6  7 ]
% [I  A2 B2 E2 mA mB mI]

function dx = moskon_simpler(x, t, args, initargs) 

    I = x(1); A = x(2); B = x(3); E = x(4);
    mA = x(5); mB = x(6); mI = x(7);

    argsCell = num2cell(args);
    [k1 k2 k3 k4 k5 k6 k_1 k_2 k_3 k_5 k_6 k_4 k7 k8 k9 k_7 k_8 k_9 b_a b_b b_i r r_a a_a a_b a_i g_a g_b g_i d_a d_b d_i d_e] = argsCell{:};

    initargsCell = num2cell(initargs);
    [e_start e_len e_d pA_init pB_init pI_init] = initargsCell{:};
  
    K1 = k1/k_1;
    K2 = k2/k_2;
    K3 = k3/k_3;
    K4 = k4/k_4;
    K5 = k5/k_5;
    K6 = k6/k_6;
    K7 = k7/k_7;
    K8 = k8/k_8;
    K9 = k9/k_9;

    enumA = (1 + r_a * K7 * A * A * K1);
    denomA = (1 + K7 * A * A * K1 + K8 * B * B * K2);
    
    enumB = (1 + r * K9 * E * E * K4);
    denomB = (1 + K9 * E * E * K4 + K7 * A * A * K3);
    
    enumI = (1 + r * K8 * B * B * K5);
    denomI = (1 + K8 * B * B * K5 + K6 * I);

    dmA = b_a * pA_init * (enumA/denomA) - a_a * mA;
    dmB = b_b * pB_init * (enumB/denomB) - a_b * mB;
    dmI = b_i * pI_init * (enumI/denomI) - a_i * mI;
    
    dA = g_a * mA - d_a * A;
    dB = g_b * mB - d_b * B;
    dI = g_i * mI - d_i * I;

    if ((t > e_start) && (t < e_start + 0.2))
      dE = e_d * 5;
    elseif ((t > e_start) && (t < e_start + e_len))
      dE = 0;
    else
      dE = - d_e * E;
    end;

    dx = [dI dA dB dE dmA dmB dmI]';
end
