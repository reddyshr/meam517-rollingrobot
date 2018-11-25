function [g] = toy_top_obj(z, N, nx, nu, dt, xgoal, ygoal)
%TRAJECTORY_COST(z) computes the cost and cost jacobian.
%   @param z: vector of decision variables containing the x_i and u_i.
%   @param N: number of sample points
%   @param nx: dimension of state vector, x
%   @param nu: dimension of input vector, u
%   @param dt: \Delta t, the inter-sample interval duration

%   @output g: total accrued cost.
%   @output dG_i: jacobian of total  accrued cost.

g = 0;
dG = zeros(N*(nx + nu),1);

% TODO: calculate g and dG
for i=1:(N-1)
   u_i_inds = (1:nu) + nx * i + nu * (i - 1);
   u_ip1_inds = (1:nu) + nx * (i+1) + nu * i;
   
   x_i_inds = (1:nx) + (nx + nu) * (i - 1);
   x_ip1_inds = x_i_inds + nx + nu;
   
   u_i = z(u_i_inds);
   u_ip1 = z(u_ip1_inds);
   
   x_i = z(x_i_inds);
   x_ip1 = z(x_ip1_inds);
   
  % g = g + (dt/2)*(dot(u_i, u_i) + dot(u_ip1, u_ip1));
  goal_vect = zeros(length(x_i(:,1)),1);
  goal_vect(4) = xgoal;
  goal_vect(5) = ygoal;
  
  x_i_diff = x_i - goal_vect;
  x_ip1_diff = x_ip1 - goal_vect;
  g = g + (dt/2)*(dot(x_i_diff, x_i_diff) + dot(x_ip1_diff, x_ip1_diff) + dot(u_i, u_i) + dot(u_ip1, u_ip1));
   
%    if (i == 1)
%        dg_ui = dt*u_i;
%        dg_uip1 = 2*dt*u_ip1;
%        
%    elseif (i == (N - 1))
%        dg_ui = 2*dt*u_i;
%        dg_uip1 = dt*u_ip1;
%        
%    else 
%        dg_ui = 2*dt*u_i;
%        dg_uip1 = 2*dt*u_ip1;
%        
%    end
%    
%    uiStart = (i-1)*(nx+nu)+1+nx;
%    uiEnd = uiStart + (nu - 1);
%    uip1Start = uiStart + (nx + nu);
%    uip1End = uiEnd + (nx + nu);
%    
%    dG(uiStart:uiEnd, 1) = dg_ui;
%    dG(uip1Start:uip1End, 1) = dg_uip1;
%        

    
end


end