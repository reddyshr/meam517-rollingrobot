
function [x_dc, u, z, N, dt, fval, INFO] = top_trajgen(M, I1, I2, I3, r, l, g, z_init, N)

   % N = 10;
    dt = 1.5/N;
    nx = 8;
    nu = 6;

    theta = pi/6;
    p = 72.84;
    s = 40;

    % Generalized Coordinates
    %q4 = x, q5 = y

    qi = zeros(5,1);
    %qi(2) = theta;

    % Generalized Speeds

    ui = zeros(3,1);
    %ui(1) = s+p*cos(theta);
    %u(3) = p*sin(theta);

    % Final COnstraints

    qf = qi;
    %qf(3) = 
    qf(4) = -0.03*2*3.14*5;
    uf = ui;


    x_0 = [qi; ui];
    x_f = [qf; uf];

    %z = find_top_trajectory(x_0, x_f, N, dt);

    [z,fval, INFO]=toy_top_min(x_0, x_f, N, dt, M, I1, I2, I3, r, l, g, z_init);

    u_t = 0:dt:(dt*(N-1));
    u = zeros(nu,0);
    x_dc = zeros(nx,0);
    for i=1:N
       x_i_inds = (1:nx) + (nx + nu) * (i - 1);
       u_i_inds = (1:nu) + nx * i + nu * (i - 1);

       x_dc(:,i) = z(x_i_inds);
       u(:,i) = z(u_i_inds);
    end

end
