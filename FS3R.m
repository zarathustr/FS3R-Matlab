% Fast Symbolic 3D Registration using MATLAB
% proposed by Jin Wu, Ming Liu, Zebo Zhou and Rui Li
%
% author: Jin Wu
% email: jin_wu_uestc@hotmail.com
% copyright 2018


function [R, T] = FS3R( b_base, r_base, weights )
    mean_b = 0;
    mean_r = 0;
    for i = 1 : length(b_base(1, :))
        mean_b = mean_b + weights(i) * b_base(:, i);
        mean_r = mean_r + weights(i) * r_base(:, i);
    end
    
    MM = 0;
    for i = 1 : length(b_base(1, :))
        MM = MM + weights(i) * (r_base(:, i) - mean_r) * (b_base(:, i) - mean_b)';
    end
    
    Hx1 = MM(1, 1);    Hx2 = MM(1, 2);    Hx3 = MM(1, 3);
    Hy1 = MM(2, 1);    Hy2 = MM(2, 2);    Hy3 = MM(2, 3);
    Hz1 = MM(3, 1);    Hz2 = MM(3, 2);    Hz3 = MM(3, 3);
    
    W = [Hx1 + Hy2 + Hz3, -Hy3 + Hz2, -Hz1 + Hx3, -Hx2 + Hy1;
        -Hy3 + Hz2, Hx1 - Hy2 - Hz3, Hx2 + Hy1, Hx3 + Hz1;
        -Hz1 + Hx3, Hx2 + Hy1, Hy2 - Hx1 - Hz3, Hy3 + Hz2;
        -Hx2 + Hy1, Hx3 + Hz1, Hy3 + Hz2, Hz3 - Hy2 - Hx1];
    
    c = det(W);
    b = 8 * Hx3 * Hy2 * Hz1 - 8 * Hx2 * Hy3 * Hz1 - ...
       8 * Hx3 * Hy1 * Hz2 + 8 * Hx1 * Hy3 * Hz2 + ...
       8 * Hx2 * Hy1 * Hz3 - 8 * Hx1 * Hy2 * Hz3;
    a = -2 * Hx1 * Hx1 - 2 * Hx2 * Hx2 - 2 * Hx3 * Hx3 - ...
        2 * Hy1 * Hy1 - 2 * Hy2 * Hy2 - 2 * Hy3 * Hy3 - ...
        2 * Hz1 * Hz1 - 2 * Hz2 * Hz2 - 2 * Hz3 * Hz3;

    T0 = 2 * a^3 + 27 * b^2 - 72 * a * c;
    theta = atan2(sqrt(4 * (a * a + 12 * c)^3 - T0 * T0), T0);
    aT1 = 1.259921049894873 * sqrt(a * a + 12 * c) * cos(theta / 3);
    T2 = sqrt( - 4 * a + 3.174802103936399 * aT1);
    lambda = 0.204124145231932 * (T2 + sqrt( - T2 * T2 - 12 * a - 29.393876913398135 * b / T2));
    
    G11 = W(1, 1) - lambda; G12 = W(1, 2); G13 = W(1, 3); G14 = W(1, 4);
    G22 = W(2, 2) - lambda; G23 = W(2, 3); G24 = W(2, 4); 
    G33 = W(3, 3) - lambda; G34 = W(3, 4); 

    q = [
		       G14 * G23 * G23 - G13 * G23 * G24 - G14 * G22 * G33 + G12 * G24 * G33 + G13 * G22 * G34 - G12 * G23 * G34;
                       G13 * G13 * G24 + G12 * G14 * G33 - G11 * G24 * G33 + G11 * G23 * G34 - G13 * G14 * G23 - G13 * G12 * G34;
	               G13 * G14 * G22 - G12 * G14 * G23 - G12 * G13 * G24 + G11 * G23 * G24 + G12 * G12 * G34 - G11 * G22 * G34;
	               - (G13 * G13 * G22 - 2 * G12 * G13 * G23 + G11 * G23 * G23 + G12 * G12 * G33 - G11 * G22 * G33)];
    q = q ./ norm(q);
    
    R = quat2dcm(q')';
    T = mean_r - R * mean_b;
end