function gen_f(mass, len, radius, grav, moi1, moi2, moi3)

    % state = [q1 q2 q3 q4 q5 u1 u2 u3]
    
    syms t q1 q2 q3 q4 q5 u1 u2 u3 tau1 tau2 tau3 f1 f2 f3 l r M g I1 I2 I3 real
  %  syms q1p1 q2p1 q3p1 q4p1 q5p1 u1p1 u2p1 u3p1 tau1p1 tau2p1 tau3p1 f1p1 f2p1 f3p1 dt real
    
    s1 = sin(q1);
    s2 = sin(q2);
    s3 = sin(q3);

    c1 = cos(q1);
    c2 = cos(q2);
    c3 = cos(q3);
    
    
    [A,b] = getTopDynamics();
    
    udot = simplify(A\b, 'Steps', 650);
    
    f = [q1; q1; q1; q1; q1; q1; q1; q1];
    
    f(1) = (c3*u1 - s3*u2)/c2;
    f(2) = (s3*u1 + c3*u2);
    f(3) = (s3*u2 - c3*u1)*(s2/c2) + u3;
    f(4) = r*(u1*(s1*s3 - c1*s2*c3) + u2*(s1*c3 + c1*s2*s3) + u3*c1*c2);
    f(5) = -r*(u1*(c1*s3 + s1*s2*c3) + u2*(c1*c3 - s1*s2*s3) - u3*s1*c2);
    f(6) = udot(1);
    f(7) = udot(2);
    f(8) = udot(3);
    
    f = subs(f, [M, l, r, g, I1, I2, I3], [mass, len, radius, grav, moi1, moi2, moi3]);
    
    matlabFunction(f, 'File', 'F', 'vars', [q1 q2 q3 q4 q5 u1 u2 u3 tau1 tau2 tau3 f1 f2 f3]);
    
%     df_dx = simplify(jacobian(f, [q1 q2 q3 q4 q5 u1 u2 u3]), 'Steps', 250);
%     
%     df_du = simplify(jacobian(f, [tau1 tau2 tau3 f1 f2 f3]), 'Steps', 250);
%     
% 
%     x_i = [q1 q2 q3 q4 q5 u1 u2 u3]';
%     x_ip1 = [q1p1 q2p1 q3p1 q4p1 q5p1 u1p1 u2p1 u3p1]';
%     u_i = [tau1 tau2 tau3 f1 f2 f3]';
%     u_ip1 = [tau1p1 tau2p1 tau3p1 f1p1 f2p1 f3p1]';
%     f_ip1 = subs(f, [q1 q2 q3 q4 q5 u1 u2 u3 tau1 tau2 tau3 f1 f2 f3], [q1p1 q2p1 q3p1 q4p1 q5p1 u1p1 u2p1 u3p1 tau1p1 tau2p1 tau3p1 f1p1 f2p1 f3p1]);
%     
%     x_iphalf = 0.5*(x_i + x_ip1) - (dt/8)*(f_ip1 - f);
%     u_iphalf = 0.5*(u_i + u_ip1);
%     
%     
%     
%     fphalf = subs(f, [q1, q2, q3, q4, q5 u1 u2 u3 tau1 tau2 tau3 f1 f2 f3], ... 
%                      [x_iphalf(1), x_iphalf(2), x_iphalf(3), x_iphalf(4), ... 
%                       x_iphalf(5), x_iphalf(6), x_iphalf(7), x_iphalf(8), ...
%                       u_iphalf(1), u_iphalf(2), u_iphalf(3), u_iphalf(4), ...
%                       u_iphalf(5), u_iphalf(6)]);
%                   
%     h_i = (3/(2*dt))*(x_ip1 - x_i) - 0.25*(f + f_ip1) - fphalf;
%     
%     dhi_dx = jacobian(h_i, [q1 q2 q3 q4 q5 u1 u2 u3]);
%     
%     dhi_dx1p = jacobian(h_i, [q1p1 q2p1 q3p1 q4p1 q5p1 u1p1 u2p1 u3p1]);
%     
%     dhi_du = jacobian(h_i, [tau1 tau2 tau3 f1 f2 f3]);
%     
%     dhi_du1p = jacobian(h_i, [q1p1 q2p1 q3p1 q4p1 q5p1 u1p1 u2p1 u3p1]);
%                   
%     dfphalf_dx = jacobian(fphalf, [q1 q2 q3 q4 q5 u1 u2 u3]);
%     dfphalf_du = jacobian(fphalf, [tau1 tau2 tau3 f1 f2 f3]);
%     
%     matlabFunction(dfphalf_dx(1,:), 'File', 'dfphalf_dx_row1', 'vars', [q1 q2 q3 q4 q5 u1 u2 u3 tau1 tau2 tau3 f1 f2 f3 ...
%         q1p1 q2p1 q3p1 q4p1 q5p1 u1p1 u2p1 u3p1 tau1p1 tau2p1 tau3p1 f1p1 f2p1 f3p1 dt]);
%     
%     matlabFunction(dfphalf_dx(2,:), 'File', 'dfphalf_dx_row2', 'vars', [q1 q2 q3 q4 q5 u1 u2 u3 tau1 tau2 tau3 f1 f2 f3 ...
%         q1p1 q2p1 q3p1 q4p1 q5p1 u1p1 u2p1 u3p1 tau1p1 tau2p1 tau3p1 f1p1 f2p1 f3p1 dt]);
%     
%     matlabFunction(dfphalf_dx(3,:), 'File', 'dfphalf_dx_row3', 'vars', [q1 q2 q3 q4 q5 u1 u2 u3 tau1 tau2 tau3 f1 f2 f3 ...
%         q1p1 q2p1 q3p1 q4p1 q5p1 u1p1 u2p1 u3p1 tau1p1 tau2p1 tau3p1 f1p1 f2p1 f3p1 dt]);
%     
%     matlabFunction(dfphalf_dx(4,:), 'File', 'dfphalf_dx_row4', 'vars', [q1 q2 q3 q4 q5 u1 u2 u3 tau1 tau2 tau3 f1 f2 f3 ...
%         q1p1 q2p1 q3p1 q4p1 q5p1 u1p1 u2p1 u3p1 tau1p1 tau2p1 tau3p1 f1p1 f2p1 f3p1 dt]);
%     
%     matlabFunction(dfphalf_dx(5,:), 'File', 'dfphalf_dx_row5', 'vars', [q1 q2 q3 q4 q5 u1 u2 u3 tau1 tau2 tau3 f1 f2 f3 ...
%         q1p1 q2p1 q3p1 q4p1 q5p1 u1p1 u2p1 u3p1 tau1p1 tau2p1 tau3p1 f1p1 f2p1 f3p1 dt]);
%     
%     matlabFunction(dfphalf_dx(6,:), 'File', 'dfphalf_dx_row6', 'vars', [q1 q2 q3 q4 q5 u1 u2 u3 tau1 tau2 tau3 f1 f2 f3 ...
%         q1p1 q2p1 q3p1 q4p1 q5p1 u1p1 u2p1 u3p1 tau1p1 tau2p1 tau3p1 f1p1 f2p1 f3p1 dt]);
%     
%      matlabFunction(dfphalf_dx(5,:), 'File', 'dfphalf_dx_row7', 'vars', [q1 q2 q3 q4 q5 u1 u2 u3 tau1 tau2 tau3 f1 f2 f3 ...
%         q1p1 q2p1 q3p1 q4p1 q5p1 u1p1 u2p1 u3p1 tau1p1 tau2p1 tau3p1 f1p1 f2p1 f3p1 dt]);
%     
%     matlabFunction(dfphalf_dx(6,:), 'File', 'dfphalf_dx_row8', 'vars', [q1 q2 q3 q4 q5 u1 u2 u3 tau1 tau2 tau3 f1 f2 f3 ...
%         q1p1 q2p1 q3p1 q4p1 q5p1 u1p1 u2p1 u3p1 tau1p1 tau2p1 tau3p1 f1p1 f2p1 f3p1 dt]);
% 
%     
%     dfphalf_dxp1 = jacobian(fphalf, [q1p1 q2p1 q3p1 q4p1 q5p1 u1p1 u2p1 u3p1]);
%     dfphalf_dup1 = jacobian(fphalf, [tau1p1 tau2p1 tau3p1 f1p1 f2p1 f3p1]);
%     
%     dfphalf_dxp1 = subs(dfphalf_dxp1, [q1p1 q2p1 q3p1 q4p1 q5p1 u1p1 u2p1 u3p1], [q1 q2 q3 q4 q5 u1 u2 u3]);
    

    

end