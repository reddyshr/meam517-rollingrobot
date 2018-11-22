function udot = udotfunc(u, q, A, b,  Mass, I_1, I_2, I_3, radius, len, grav)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    syms t q1 q2 q3 u1 u2 u3 u1_dot u2_dot u3_dot x_dot y_dot  l r M g tau f I1 I2 I3 real
    
    p = double(subs(A, [q1 q2 q3 u1 u2 u3 r l M g tau f I1 I2 I3], [q(1) q(2) q(3) u(1) u(2) u(3) radius len Mass grav tau f I_1 I_2 I_3]));
    n = double(subs(b, [q1 q2 q3 u1 u2 u3 r l M g tau f I1 I2 I3], [q(1) q(2) q(3) u(1) u(2) u(3) radius len Mass grav tau f I_1 I_2 I_3]));

    udot = p \ n;


end