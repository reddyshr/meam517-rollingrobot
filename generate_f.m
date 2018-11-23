function generate_f()

    syms t q1 q2 q3 u1 u2 u3 q1_t(t) q2_t(t) q3_t(t) u1_t(t) u2_t(t) u3_t(t) q1_f q2_f q3_f q1_dot_f q2_dot_f q3_dot_f q1_ddot_f q2_ddot_f q3_ddot_f l r M g tau f I1 I2 I3 real
    
    s1_t = sin(q1_t(t));
    s2_t = sin(q2_t(t));
    s3_t = sin(q3_t(t));

    c1_t = cos(q1_t(t));
    c2_t = cos(q2_t(t));
    c3_t = cos(q3_t(t));
    
    q1_dot = diff(q1_t(t), t);
    q2_dot = diff(q2_t(t), t);
    q3_dot = diff(q3_t(t), t);
    
    u1_exp = c2_t*c3_t*q1_dot + ((1 - c3_t^2) / s3_t)*q2_dot;
    u2_exp = c3_t*q2_dot - q1_dot*c2_t*c3_t;
    u3_exp = q3_dot - (s3_t*u2_exp - c3_t*u1_exp) * (s2_t/c2_t);
    
    u1_dot_exp = diff(u1_exp, t);
    u2_dot_exp = diff(u2_exp, t);
    u3_dot_exp = diff(u3_exp, t);
    u_dot = [q1; q1; q1];
    u_dot(1) = u1_dot_exp;
    u_dot(2) = u2_dot_exp;
    u_dot(3) = u3_dot_exp;
    
    [A,b] = getTopDynamics();
    
    A = subs(A, [q1 q2 q3 u1 u2 u3], [q1_t(t) q2_t(t) q3_t(t) u1_exp u2_exp u3_exp]);
    b = subs(b, [q1 q2 q3 u1 u2 u3], [q1_t(t) q2_t(t) q3_t(t) u1_exp u2_exp u3_exp]);
    
    eom = A*u_dot - b;
    eom = subs(eom, [q1_t(t) q2_t(t) q3_t(t) diff(q1_t(t), t) diff(q2_t(t), t) diff(q3_t(t), t) diff(diff(q1_t(t), t), t) diff(diff(q2_t(t), t), t) diff(diff(q3_t(t), t), t)],  [q1_f q2_f q3_f q1_dot_f q2_dot_f q3_dot_f q1_ddot_f q2_ddot_f q3_ddot_f]);
    
    eqns = [eom(1) == 0; eom(2) == 0; eom(3) == 0];
    
    %eom = simplify(eom == [0; 0; 0],  'Steps', 500);
    
    [q1ddot_exp, q2ddot_exp, q3ddot_exp] = solve(eqns, [q1_ddot_f; q2_ddot_f; q3_ddot_f]);

   j = 2;

end