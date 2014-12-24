% x =
% [1  2  3  4  5  6  7  8  9  10]
% [A  B  E  I  A2 B2 E2 mA mB mI]

function dx = moskon_simple(x, t, args) 

    A = x(1); B = x(2); E = x(3); I = x(4);
    A2 = x(5); B2 = x(6); E2 = x(7);
    mA = x(8); mB = x(9); mI = x(10);
    
    argsCell = num2cell(args);
    [e_switch_point e_switch_d pA_init pB_init pI_init k1 k2 k3 k4 k5 k6 k_1 k_2 k_3 k_5 k_6 k_4 k7 k8 k9 k_7 k_8 k_9 b_a b_b b_i fi fi_a a_a a_b a_i g_a g_b g_i d_a d_b d_i d_e] = argsCell{:};

    % Fractional occupancy for pA %
    w0 = 1;
    w1 = A2/(k_1/k1);
    w2 = B2/(k_2/k2);
    wsum = w0+w1+w2;
    
    p0 = w0/wsum;
    p1 = w1/wsum;
    p2 = w2/wsum;

    pA   = p0*pA_init;
    pAA2 = p1*pA_init;
    pAB2 = p2*pA_init;
    
    % Fractional occupancy for pB %
    w1 = A2/(k_3/k3);
    w2 = E2/(k_4/k4);
    wsum = w0+w1+w2;
    
    p0 = w0/wsum;
    p1 = w1/wsum;
    p2 = w2/wsum;

    pB   = p0*pB_init;
    pBA2 = p1*pB_init;
    pBE2 = p2*pB_init;
	    
    % Fractional occupancy for pI %
    w1 = B2/(k_5/k5);
    w2 = I/(k_6/k6);
    wsum = w0+w1+w2;
    
    p0 = w0/wsum;
    p1 = w1/wsum;
    p2 = w2/wsum;

    pI   = p0*pI_init;
    pIB2 = p1*pI_init;
    pII  = p2*pI_init;

    % the rest of equations %
    dx(1) = -2*k7*A^2 + 2*k_7*A2 + g_a*mA - d_a*A;
    dx(2) = -2*k8*B^2 + 2*k_8*B2 + g_b*mB - d_b*B;

    %switch tabs with spaces!
    if (t > e_switch_point) && (t < e_switch_point+1) 
      dx(3) = e_switch_d;
    else
      dx(3) = -2*k9*E^2 + 2*k_9*E2 - d_e*E;
    end

    dx(4) = -k6*I*pI  + k_6*pII  + g_i*mI - d_i*I;
    
    dx(5) = -k1*pA*A2 + k_1*pAA2 - k3*A2*pB  + k_3*pBA2 + k7*A^2 - k_7*A2 - d_a*A2;
    dx(6) = -k2*pA*B2 + k_2*pAB2 - k5*B2*pI  + k_5*pIB2 + k8*B^2 - k_8*B2 - d_b*B2;
    dx(7) = -k4*pB*E2 + k_4*pBE2                        + k9*E^2 - k_9*E2 - d_e*E2;
    
    dx(8)  = b_a*pA + b_a*fi_a*pAA2 - a_a*mA;
    dx(9)  = b_b*pB + b_b*fi  *pBE2 - a_b*mB;
    dx(10) = b_i*pI + b_i*fi  *pIB2 - a_i*mI;
end
