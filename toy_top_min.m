function [z,fval, INFO]=toy_top_min(x_0, x_f, N, dt, M, I1, I2, I3, r, l, g, z_init)
% Mimics sntoyA.f in $SNOPT/examples
% Example of the 'fmincon'-style call to SNOPT.
%
%     Minimize      3*x(1) + (x(1) + x(2) + x(3))^2 + 5*x(4)
%
%     subject to             4*x(2)   + 2*x(3)               >= 0
%                     x(1) +   x(2)^2 +   x(3)^2              = 2
%                              x(2)^4 +   x(3)^4   +   x(4)   = 4
%
%                     x(1) >= 0,                       x(4) >= 0.
%
%

snscreen on;
snprint('toymin.out');  % By default, screen output is off;

sntoy.spc = which('sntoy.spc');
snspec (sntoy.spc);

snseti ('Major Iteration limit', 250);

  nx = 8;
  nu = 6;
  
  xgoal = -1;
  ygoal = 1;
  
  % TODO: Add constraints to Aeq, beq to enforce starting at x_0 and ending
  % at x_f
  x_0_inds = 1:nx;
  x_f_inds = x_0_inds + (N - 1) * (nx + nu);
  Aeq = zeros(2*nx, N * (nx + nu));
  beq = zeros(2*nx, 1);
  
  beq = [x_0; 0; 0; 0; xgoal; ygoal; 0; 0; 0];
 % beq = [x_0; x_f];
  Aeq(1:nx, 1:nx) = eye(nx);
   eq_conds = zeros(nx,nx);
%  % eq_conds(1:3,1:3) = eye(3);
   eq_conds(4,4) = 1;
   eq_conds(5,5) = 1;
   eq_conds(6,6) = 1;
   eq_conds(7,7) = 1;
   eq_conds(8,8) = 1;
 %  eq_cond(8,8) = 1;
   Aeq(nx+1:end, (end-nu-nx+1):(end-nu)) = eq_conds;
  
  A = zeros(2*nx, N * (nx + nu));
  b = zeros(2*nx, 1);
  cond = zeros(nx*2, nx);
%   cond(1,4) = -1;
%   cond(2,4) = 1;
%   cond(3,5) = -1;
%   cond(4,5) = 1;
%   cond(5,6) = -1;
%   cond(6,6) = 1;
%   cond(7,7) = -1;
%   cond(8,7) = 1;
%   cond(9,8) = -1;
%   cond(10,8) = 1;
  
%   cond = [0 0 0 -1 0 0 0 0;
%           0 0 0 1 0 0 0 0;
%           0 0 0 0 -1 0 0 0;
%           0 0 0 0 1 0 0 0;
%           0 0 0 0 0 -1 0 0;
%           0 0 0 0 0 1 0 0;
%           0 0 0 0 0 0 -1 0;
%           0 0 0 0 0 0 1 0];
  b = [0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0; 0];
  A(1:end, (end-nu-nx+1):(end-nu)) = cond;
  
  
  MaxU = 50;
  
  % TODO: Add bounding box constraints u_1 \in [-M,M], u_2 \in [-M,M]
  lb = -inf(N * (nx + nu),1);
  ub = inf(N * (nx + nu),1);
  for i=1:N
      u_i_inds = (1:nu) + nx * i + nu * (i - 1);
      lb(u_i_inds) = -MaxU;
      ub(u_i_inds) = MaxU;
      
%       lb(u_i_inds(1)) = 0;
%       ub(u_i_inds(1)) = 0;
%       lb(u_i_inds(2)) = 0;
%       ub(u_i_inds(2)) = 0;
      
      lb(u_i_inds(4)) = 0;
      ub(u_i_inds(4)) = 0;
      lb(u_i_inds(5)) = 0;
      ub(u_i_inds(5)) = 0;
      lb(u_i_inds(6)) = 0;
      ub(u_i_inds(6)) = 0;


  end
  
  
  % TODO: make initial guess for z
  z0 = zeros(N * (nx + nu), 1);
  
  if (z_init == 0)
      slope = (x_f - x_0)/N;
      for i=1:N
          u_i_inds = (1:nu) + nx * i + nu * (i - 1);
          x_i_inds = (1:nx) + (nx + nu) * (i - 1);
          z0(x_i_inds) = x_0 + (i - 1)*slope;
          z0(u_i_inds(3)) = -0.1;
          z0(u_i_inds(2)) = 0.1;
      end
  else
      z0 = z_init;
  end
  
options.name = 'topprob';
options.stop = @toySTOP;

[z,fval,INFO,output,lambda,states] = snsolve( @(z) toy_top_obj(z, N, nx, nu, dt, xgoal, ygoal), z0, A, b, Aeq, beq, lb, ub, ...
					      @(z) toy_top_con(z, N, nx, nu, dt, M, I1, I2, I3, r, l, g), options);

snprint off;
snend;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [iAbort] = toySTOP(itn, nMajor, nMinor, condZHZ, obj, merit, step, ...
			      primalInf, dualInf, maxViol, maxViolRel, ...
			      x, xlow, xupp, xmul, xstate, ...
			      F, Flow, Fupp, Fmul, Fstate)

% Called every major iteration
% Use iAbort to stop SNOPT (if iAbort == 0, continue; else stop)

iAbort = 0