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
    
    W = [MM(1, 1) + MM(2, 2) + MM(3, 3), -MM(2, 3) + MM(3, 2), -MM(3, 1) + MM(1, 3), -MM(1, 2) + MM(2, 1);
        -MM(2, 3) + MM(3, 2), MM(1, 1) - MM(2, 2) - MM(3, 3), MM(1, 2) + MM(2, 1), MM(1, 3) + MM(3, 1);
        -MM(3, 1) + MM(1, 3), MM(1, 2) + MM(2, 1), MM(2, 2) - MM(1, 1) - MM(3, 3), MM(2, 3) + MM(3, 2);
        -MM(1, 2) + MM(2, 1), MM(1, 3) + MM(3, 1), MM(2, 3) + MM(3, 2), MM(3, 3) - MM(2, 2) - MM(1, 1)];
    
    c = det(W);
    b = 8 * MM(1, 3) * MM(2, 2) * MM(3, 1) - 8 * MM(1, 2) * MM(2, 3) * MM(3, 1) - ...
       8 * MM(1, 3) * MM(2, 1) * MM(3, 2) + 8 * MM(1, 1) * MM(2, 3) * MM(3, 2) + ...
       8 * MM(1, 2) * MM(2, 1) * MM(3, 3) - 8 * MM(1, 1) * MM(2, 2) * MM(3, 3);
    a = -2 * MM(1, 1) * MM(1, 1) - 2 * MM(1, 2) * MM(1, 2) - 2 * MM(1, 3) * MM(1, 3) - ...
        2 * MM(2, 1) * MM(2, 1) - 2 * MM(2, 2) * MM(2, 2) - 2 * MM(2, 3) * MM(2, 3) - ...
        2 * MM(3, 1) * MM(3, 1) - 2 * MM(3, 2) * MM(3, 2) - 2 * MM(3, 3) * MM(3, 3);

    T0 = 2 * a * a * a + 27 * b * b - 72 * a * c;
    theta = atan2(sqrt(4 * (a * a + 12 * c)^3 - T0 * T0), T0);
    aT1 = sqrt(a * a + 12 * c) * cos(theta / 3);
    T2 = 2 * sqrt(aT1 - a);
    
    if(abs(b) < 1e-5 && abs(T2) < 1e-5)
        lambda = sqrt( - b / 2);
    else
        lambda = 0.204124145231932 * (T2 + sqrt( - T2 * T2 - 12 * a - 29.393876913398135 * b / T2));
    end
    
    G11 = W(1, 1) - lambda; G12 = W(1, 2); G13 = W(1, 3); G14 = W(1, 4);
    G22 = W(2, 2) - lambda; G23 = W(2, 3); G24 = W(2, 4); 
    G33 = W(3, 3) - lambda; G34 = W(3, 4); 

    q = [
		       G14 * G23 * G23 - G13 * G23 * G24 - G14 * G22 * G33 + G12 * G24 * G33 + G13 * G22 * G34 - G12 * G23 * G34;
                       G13 * G13 * G24 + G12 * G14 * G33 - G11 * G24 * G33 + G11 * G23 * G34 - G13 * G14 * G23 - G13 * G12 * G34;
	               G13 * G14 * G22 - G12 * G14 * G23 - G12 * G13 * G24 + G11 * G23 * G24 + G12 * G12 * G34 - G11 * G22 * G34;
	               - (G13 * G13 * G22 - 2 * G12 * G13 * G23 + G11 * G23 * G23 + G12 * G12 * G33 - G11 * G22 * G33)];
    R = quat2dcm(q')';
    T = mean_r - R * mean_b;
end
