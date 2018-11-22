function runsim()
    close all;

    M = 1;
    I1 = 1.25*10^-3;
    I2 = 6.25*10^-4;
    I3 = 6.25*10^-4;
    r = 0.02;
    l = 0.03;
    g = 9.81;

    theta = pi/6;
    p = 72.84;
    s = 50;

    % Generalized Coordinates
    %q4 = x, q5 = y

    q = zeros(5,1);
    q(2) = theta;

    % Generalized Speeds

    u = zeros(5,1);
    u(1) = s+p*cos(theta);
    u(3) = p*sin(theta);
    
    % Simulation Parameters

    TMAX = 30;
    STEP = 0.001;
    time = 0;
    cg3arr = [];

    figure('units','normalized','outerposition',[0 0 1 1]);

    while (time < TMAX)
        u(4) = 0;
        u(5) = 0;

        tspan = [time time+STEP];

        [~, qtemp] = ode45(@(time,y) qdotfunc(y, u, r), tspan, q);
        
        
        f = @(q, u) A(q(1),q(2),q(3),I1,I2,I3,M,l,r) \ b(q(1),q(2),q(3),u(1),u(2),u(3),I1,I2,I3,M,l,r,g);
        
        [~, utemp] = ode45(@(t, y) f(q, y), tspan, u(1:3));
        
        q = wrapToPi(qtemp(end, :))'; 
        u(1:3) = utemp(end,:)';

        s1 = sin(q(1));
        s2 = sin(q(2));
        s3 = sin(q(3));

        c1 = cos(q(1));
        c2 = cos(q(2));
        c3 = cos(q(3));


        H = [c1*c3-s1*s2*s3, -c2*s1, c1*s3+c3*s1*s2;
             c3*s1+c1*s2*s3, c1*c2, s1*s3-c1*c3*s2;
             -c2*s3,         s2,    c2*c3];

        bz = H * [0; 0; 0.01];
        bx = H * [0.01; 0; 0];


        bz_unit = bz ./ norm(bz);
        cg = (l+r)*bz_unit;
        cg(1) = cg(1) + q(4);
        cg(2) = cg(2) + q(5);
        cg3arr = [cg3arr cg(3)];
        p_ag = [cg(1) - q(4); cg(2) - q(5); cg(3)];




        plot3([cg(1) cg(1)+bz(1)], [cg(2) cg(2)+bz(2)], [cg(3) cg(3)+bz(3)], 'k');
        hold on;
        xlim([-1 1]);
        ylim([-1 1]);
        zlim([0 0.2]);
        grid ON;
        plot3([cg(1) cg(1)+bx(1)], [cg(2) cg(2)+bx(2)], [cg(3) cg(3)+bx(3)], 'r');
        plot3([q(4) q(4)+p_ag(1)],[q(5) q(5)+p_ag(2)],[0 p_ag(3)], '-m');
        hold off;
        pause(0.0005);

        time = time + STEP;

        disp(time);

    end
end