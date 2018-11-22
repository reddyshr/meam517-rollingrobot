function [A,b] = getTopDynamics()
    
% Returns the dynamics of the system in the form of Ax = b, where x is
% a vector of the time derivatives of the generalized speeds. Since we are
% only considering rolling, x = [u1_dot u2_dot u3_dot]

    syms t q1 q2 q3 u1 u2 u3 u1_dot u2_dot u3_dot x_dot y_dot  l r M g tau f I1 I2 I3 real

    %% System Setup

    s1 = sin(q1);
    s2 = sin(q2);
    s3 = sin(q3);

    c1 = cos(q1);
    c2 = cos(q2);
    c3 = cos(q3);

    T = [c2*c3, -c2*s3, s2;
         c1*s3+c3*s1*s2, c1*c3-s1*s2*s3, -c2*s1;
         s1*s3-c1*c3*s2, c3*s1+c1*s2*s3, c1*c2];
     
    x_dot = [0; r*(u1*(s1*s3 - c1*s2*c3) + u2*(s1*c3 + c1*s2*s3) + u3*c1*c2); 0];
    y_dot = [0; 0; -r*(u1*(c1*s3 + s1*s2*c3) + u2*(c1*c3 - s1*s2*s3) - u3*s1*c2)];

    %% Kinematics of System

    w_ab = [u1; u2; u3];
    alpha_ab = [u1_dot; u2_dot; u3_dot];
    r_pg = [l; 0; 0];

    x_dot_b = T'*x_dot;
    y_dot_b = T'*y_dot;

    A_Vpg_b = cross(w_ab, r_pg);
    A_Vop_b = x_dot_b + y_dot_b;

    A_Vg_b = A_Vpg_b + A_Vop_b;

    %% Generalized Active Forces

    %Tau_b = [0; tau; 0];
    %R_b = T'*[-M*g; 0; 0] + [0; f; 0];
    r_pg_a = T*r_pg;
    momentarm = sqrt(r_pg_a(2)^2 + r_pg_a(3)^2);
    tau_vect_a = cross([1;0;0], r_pg_a);
    unit_tau_vect_a = tau_vect_a / norm(tau_vect_a);
    Tau_b = T'*((M*g*momentarm)*unit_tau_vect_a);

    R_b = [0; 0; 0];

    %Tau_b = [0; 0; 0];
    %R_b = T'*[-M*g; 0; 0];

    F1_b = dot(diff(w_ab, u1), Tau_b) + dot(diff(A_Vg_b, u1), R_b);
    F2_b = dot(diff(w_ab, u2), Tau_b) + dot(diff(A_Vg_b, u2), R_b);
    F3_b = dot(diff(w_ab, u3), Tau_b) + dot(diff(A_Vg_b, u3), R_b);

    %% Generalized Inertia Forces
% 
    Tau_inertia_b = [-(alpha_ab(1)*I1 - u2*u3*(I2 - I3));
                     -(alpha_ab(2)*I2 - u3*u1*(I3 - I1));
                     -(alpha_ab(3)*I3 - u1*u2*(I1 - I2))];


    A_Ag_a = [u1; u1; u1];
    A_Ag_a(1) = l*(-u2_dot*s2-u3_dot*c2*s3-(u2^2 + ... 
        + u3^2)*c2*c3-u1*u2*c2*s3 + u3*u1*s2);

    A_Ag_a(2) = r*u1_dot*(s1*s3-s2*c1*c3) + u2_dot*(r*(s2*c1*s3 + s1*c3) + l*s1*c2) + ... 
        + u3_dot*(r*c1*c2 + l*(c1*c3-s1*s2*s3)) + ... 
        - (u2^2 + u3^2)*l*(c1*s3+s1*s2*c3) + ...
        + u1*u2*l*(c1*c3 - s1*s2*s3) - u3*u1*l*s1*c2;

    A_Ag_a(3) = -r*u1_dot*(s1*s2*c3 + c1*s3) + u2_dot*(r*(s1*s2*s3 - c1*c3)-l*c1*c2) + ...
        + u3_dot*(r*c2*s1+l*(s1*c3+c1*s2*s3)) + ...
        - (u2^2 + u3^2)*l*(s1*s3 - c1*s2*c3) + ...
        + u1*u2*l*(s1*c3+c1*s2*s3)+u3*u1*l*c1*c2;

    A_Ag_b = T'*A_Ag_a;

    R_inertia_b = -M*A_Ag_b;

    F1_inertia_b = dot(diff(w_ab, u1), Tau_inertia_b) + dot(diff(A_Vg_b, u1), R_inertia_b);
    F2_inertia_b = dot(diff(w_ab, u2), Tau_inertia_b) + dot(diff(A_Vg_b, u2), R_inertia_b);
    F3_inertia_b = dot(diff(w_ab, u3), Tau_inertia_b) + dot(diff(A_Vg_b, u3), R_inertia_b);

    eom1 = simplify(simplify(F1_b, 'Steps', 150) + simplify(F1_inertia_b,  'Steps', 150),  'Steps', 500);
    eom2 = simplify(simplify(F2_b, 'Steps', 150) + simplify(F2_inertia_b,  'Steps', 150),  'Steps', 500);
    eom3 = simplify(simplify(F3_b, 'Steps', 150) + simplify(F3_inertia_b,  'Steps', 150),  'Steps', 500);

    eqns = [eom1 == 0,
            eom2 == 0,
            eom3 == 0];
    vars = [u1_dot, u2_dot, u3_dot];

    [A, b] = equationsToMatrix(eqns, vars);
    matlabFunction(A, 'File', 'A', 'vars', [q1 q2 q3 I1 I2 I3 M l r]);
    matlabFunction(b, 'File', 'b', 'vars', [q1 q2 q3 u1 u2 u3 I1 I2 I3 M l r g]);
end

