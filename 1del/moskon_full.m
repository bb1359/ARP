% x =
% [1  2  3  4    5    6    7    8    9   ...]
% [pA pB pI pAA2 pAB2 pBA2 pBE2 pIB2 pII ...]
% [... 10 11 12 13 14 15 16 17 18 19]
% [... A  B  E  I  A2 B2 E2 mA mB mI]

function dx = moskon_full(x, t, args)

    pA = x(1); pB = x(2); pI = x(3);
    pAA2 = x(4); pAB2 = x(5);
    pBA2 = x(6); pBE2 = x(7);
    pIB2 = x(8); pII = x(9);
    A = x(10); B = x(11); E = x(12); I = x(13);
    A2 = x(14); B2 = x(15); E2 = x(16);
    mA = x(17); mB = x(18); mI = x(19);
    
    argsCell = num2cell(args);
    [e_switch_point e_switch_d pA_init pB_init pI_init k1 k2 k3 k4 k5 k6 k_1 k_2 k_3 k_5 k_6 k_4 k7 k8 k9 k_7 k_8 k_9 b_a b_b b_i fi fi_a a_a a_b a_i g_a g_b g_i d_a d_b d_i d_e] = argsCell{:};

    dx(1) = -k1*pA*A2 + k_1*pAA2 - k2*pA*B2 + k_2*pAB2;
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

    if (t > e_switch_point) && (t < e_switch_point+1)
      dx(12) = e_switch_d;
    else
      dx(12) = -2*k9*E^2 + 2*k_9*E2 - d_e*E;
    end

    dx(13) = -k6*I*pI  + k_6*pII  + g_i*mI - d_i*I;
    
    dx(14) = -k1*pA*A2 + k_1*pAA2 - k3*A2*pB  + k_3*pBA2 + k7*A^2 - k_7*A2 - d_a*A2;
    dx(15) = -k2*pA*B2 + k_2*pAB2 - k5*B2*pI  + k_5*pIB2 + k8*B^2 - k_8*B2 - d_b*B2;
    dx(16) = -k4*pB*E2 + k_4*pBE2                        + k9*E^2 - k_9*E2 - d_e*E2;
    
    dx(17) = b_a*pA + b_a*fi_a*pAA2 - a_a*mA;
    dx(18) = b_b*pB + b_b*fi  *pBE2 - a_b*mB;
    dx(19) = b_i*pI + b_i*fi  *pIB2 - a_i*mI;
end
