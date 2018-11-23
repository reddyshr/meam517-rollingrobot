function [c, ceq, dC, dCeq] = all_constraints(z, N, nx, nu, dt)

[ceq, dCeq] = dynamics_constraints(z, N, nx, nu, dt);

c = zeros(0,1);
dC = zeros(0,numel(z));

dC = sparse(dC)';
dCeq = sparse(dCeq)';
end