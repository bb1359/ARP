% x =
% [1  2  3  4  5  6  7  8  9  10]
% [A  B  E  I  A2 B2 E2 mA mB mI]

function dx = moskon_full(x, t)

    A = x(1); B = x(2); E = x(3); I = x(4);
    A2 = x(5); B2 = x(6); E2 = x(7);
    mA = x(8); mB = x(9); mI = x(10);
    
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

    e_switch_point = 25;
    e_switch_d = 1000;

    % Fractional occupancy for pA %
    pA_0 = 1;
    w0 = 1;
    w1 = A2/(k_1/k1);
    w2 = B2/(k_2/k2);
    wsum = w0+w1+w2;
    
    p0 = w0/wsum;
    p1 = w1/wsum;
    p2 = w2/wsum;

    pA = p0*pA_0;
    pAA2 = p1*pA_0;
    pAB2 = p2*pA_0;
    
    % Fractional occupancy for pB %
    pB_0 = 0.5;
    w0 = 1;
    w1 = A2/(k_3/k3);
    w2 = E2/(k_4/k4);
    wsum = w0+w1+w2;
    
    p0 = w0/wsum;
    p1 = w1/wsum;
    p2 = w2/wsum;

    pB = p0*pB_0;
    pBA2 = p1*pB_0;
    pBE2 = p2*pB_0;
	    
    % Fractional occupancy for pI %
    pI_0 = 1;
    w0 = 1;
    w1 = B2/(k_5/k5);
    w2 = I/(k_6/k6);
    wsum = w0+w1+w2;
    
    p0 = w0/wsum;
    p1 = w1/wsum;
    p2 = w2/wsum;

    pI = p0*pI_0;
    pIB2 = p1*pI_0;
    pII = p2*pI_0;

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
