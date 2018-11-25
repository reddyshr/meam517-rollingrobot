function [c,ceq] = toy_top_con(z, N, nx, nu, dt, M, I1, I2, I3, r, l, g)

ceq = 0;

[ceq, dCeq] = dynamics_constraints(z, N, nx, nu, dt, M, I1, I2, I3, r, l, g);

c = zeros(0,1);
dC = zeros(0,numel(z));

dC = sparse(dC)';
dCeq = sparse(dCeq)';

end