function [z, Aeq, beq, lb, ub, z0] = find_top_trajectory(x_0, x_f, N, dt)
%FIND_SWINGUP_TRAJECTORY(x_0, x_f, N, dt) executes a direct collocation
%optimization program to find an input sequence to drive the cartpole
%system from x_0 to x_f.
%
%   @param x_0: the state at the start of the trajectory
%   @param x_f: the state at the emd of the trajectory
%   @param N: number of state/input samples
%   @param dt: \Delta t, the duration of the interval between samples
%
%   @output z: decision variable vector containing the x_i and u_i
%   @output Aeq: matrix from linear constrant Aeq z = beq
%   @output beq: vector from linear constrant Aeq z = beq
%   @output lb: lower bound from box constraint lb <= z <= ub
%   @output ub: upper bound from box constraint lb <= z <= ub
%   @output z0: initial guess for z

  nx = 8;
  nu = 6;
  
  % TODO: Add constraints to Aeq, beq to enforce starting at x_0 and ending
  % at x_f
  x_0_inds = 1:nx;
  x_f_inds = x_0_inds + (N - 1) * (nx + nu);
  Aeq = zeros(2*nx, N * (nx + nu));
  beq = zeros(2*nx, 1);
  
 % beq = [x_0; x_f(1:3); 0; 0; 0; 0; 0];
  beq = [x_0; x_f];
  Aeq(1:nx, 1:nx) = eye(nx);
  eq_conds = zeros(nx,nx);
  eq_conds(1:3,1:3) = eye(3);
  Aeq(nx+1:end, (end-nu-nx+1):(end-nu)) = eye(nx);
  
  A = zeros(2*nx, N * (nx + nu));
  b = zeros(2*nx, 1);
  cond = [0 0 0 1 0 0 0 0;
          0 0 0 -1 0 0 0 0;
          0 0 0 0 1 0 0 0;
          0 0 0 0 -1 0 0 0;
          0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0;
          0 0 0 0 0 0 0 0];
  b = [0; 0; 0; 0; 0; 0; 0; 0; 1.1; -0.9; 0.1; -0.1; 0; 0; 0; 0];
  A(nx+1:end, (end-nu-nx+1):(end-nu)) = cond;
  
  
  M = 50;
  
  % TODO: Add bounding box constraints u_1 \in [-M,M], u_2 \in [-M,M]
  lb = -inf(N * (nx + nu),1);
  ub = inf(N * (nx + nu),1);
  for i=1:N
      u_i_inds = (1:nu) + nx * i + nu * (i - 1);
      lb(u_i_inds) = -M;
      ub(u_i_inds) = M;
      
      lb(u_i_inds(1)) = 0;
      ub(u_i_inds(1)) = 0;
      
      lb(u_i_inds(4:end)) = 0;
      ub(u_i_inds(4:end)) = 0;

  end
  
  % TODO: make initial guess for z
  z0 = zeros(N * (nx + nu), 1);
  slope = (x_f - x_0)/N;
  for i=1:N
      u_i_inds = (1:nu) + nx * i + nu * (i - 1);
      x_i_inds = (1:nx) + (nx + nu) * (i - 1);
      z0(x_i_inds) = x_0 + (i - 1)*slope;
      z0(u_i_inds(3)) = -0.1;

  end
  
  options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true, 'CheckGradients', false, 'OptimalityTolerance', 3e0, 'MaxIterations', 6000, 'MaxFunctionEvaluations', 20000, 'Display','iter');
  problem.objective = @(z) trajectory_cost(z, N, nx, nu, dt);
  
  
  problem.x0 = z0;
  problem.options = options;
  problem.nonlcon = @(z) all_constraints(z, N, nx, nu, dt);
  problem.solver = 'fmincon';
  problem.Aeq = Aeq;
  problem.beq = beq;
  %problem.A = A;
 % problem.b = b;
  problem.lb = lb;
  problem.ub = ub;

  z = fmincon(problem);
end