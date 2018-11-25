function h_i = dynamics_constraint(x_i, u_i, x_ip1, u_ip1, dt, M, I1, I2, I3, r, l, g)
%DYNAMICS_CONSTRAINT(xm, um, xp, up, dt) computes the vector contstraint
%   h_i = \dot{s}_i(\Delta t/2) - f(s_i(\Delta t/2), 0.5(u_{i+1}+u{i})) = 0
%
%   @param x_i: x_i, the state at the start of the interval.
%   @param u_i: u_i, the input at the start of the interval.
%   @param x_ip1: x_{i+1}, the state at the end of the interval.
%   @param u_ip1: u_{i+1}, the input at the end of the interval.
%   @param dt: \Delta t, the duration of the interval
%
%   @output h_i: quantity from above expresion that should be 0

% x = [q1 q2 q3 q4 q5 u1 u2 u3];
% u = [tau1 tau tau3 f1 f2 f3];

f_i = F(x_i(1), x_i(2), x_i(3), x_i(4), x_i(5), x_i(6), x_i(7), x_i(8), ...
        u_i(1), u_i(2), u_i(3), u_i(4), u_i(5), u_i(6), M,l,r,g,I1,I2,I3);
f_ip1 = F(x_ip1(1), x_ip1(2), x_ip1(3), x_ip1(4), x_ip1(5), x_ip1(6), x_ip1(7), x_ip1(8), ... 
          u_ip1(1), u_ip1(2), u_ip1(3), u_ip1(4), u_ip1(5), u_ip1(6), M,l,r,g,I1,I2,I3);

x_iphalf = 0.5*(x_i + x_ip1) - (dt/8)*(f_ip1 - f_i);
u_iphalf = 0.5*(u_i + u_ip1);
xdot_iphalf = (3/(2*dt))*(x_ip1 - x_i) - (1/4)*(f_i + f_ip1);

% TODO: compute h_i.
h_i = xdot_iphalf - F(x_iphalf(1), x_iphalf(2), x_iphalf(3), x_iphalf(4), ...
                      x_iphalf(5), x_iphalf(6), x_iphalf(7), x_iphalf(8), ...
                      u_iphalf(1), u_iphalf(2), u_iphalf(3), u_iphalf(4), ...
                      u_iphalf(5), u_iphalf(6), M,l,r,g,I1,I2,I3);

end