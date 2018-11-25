function runsim()
    close all;
    M = 1;
    I1 = 1.25*10^-3;
    I2 = 6.25*10^-4;
    I3 = 6.25*10^-4;
    r = 0.1;
    l = 0.0;
    g = 9.81;
    
    z_init = 0;
    numOfRuns = 0;
    goAgain = true;
    N = 10;
    goalN = 40;
    
    noise = load('train');
    
    while (goAgain) 
         [x_col, u_col, z_col, N, dt, fval, INFO] = top_trajgen(M, I1, I2, I3, r, l, g, z_init, N);
         z_init = z_col;
         numOfRuns = numOfRuns + 1;
         disp(u_col);
         disp(x_col);
         
         goodResponse = false;
         while (~goodResponse)             
             prompt = 'Finish(1) -- Reseed w/ same N(2) -- Reseed w/ N+10(3)?';
             response = 0;
             %response = input(prompt);
             if (response == 2 || INFO == 32)
                 goAgain = true;
                 goodResponse = true;
             elseif (response == 1 || (INFO == 1 && N == goalN)) 
                 goAgain = false;
                 goodResponse = true;
                 
             elseif (response == 3 || (INFO == 1 && N < goalN))
                 goAgain = true;
                 goodResponse = true;
                 ti = linspace(1,N,N);
                 tf = linspace(1,N,N+10);                 
                 N = N + 10;
                 x_temp = interp1(ti, x_col', tf)';
                 u_temp = interp1(ti, u_col', tf)';
                 z_temp = [];
                 for i = 1:N
                     z_temp = [z_temp; x_temp(:,i); u_temp(:,i)];
                 end
                 z_init = z_temp;
             else
                 goodResponse = true;
                 disp('Try Again');
             end
         end
         
    end
    
    sound(noise.y, noise.Fs);
    
   %%
       
    % Simulation Parameters
    
    q = zeros(5,1);
    u = zeros(3,1);
    
    q = x_col(1:5,1);
    u = x_col(6:8,1);

    TMAX = 5;
    STEP = 0.001;
    time = 0;
    cg3arr = [];
    
    control_des = getContinuousInputTrajectory(u_col, N, dt, STEP);
    control_len = length(control_des);

    save('/home/reddyshr/Desktop/Meam517/meam517-finalproject/trials/trial11.mat')
    
    %%
    close all;
    figure('units','normalized','outerposition',[0 0 1 1]);

    ind = 1;
    
    
    
    while (time < TMAX)
        u(4) = 0;
        u(5) = 0;
        
        if (ind <= control_len)
            tau1 = control_des(1,ind);
            tau2 = control_des(2,ind);
            tau3 = control_des(3,ind);
            f1 = control_des(4,ind);
            f2 = control_des(5,ind);
            f3 = control_des(6,ind);
        else 
            tau1 = 0;
            tau2 = 0;
            tau3 = 0;
            f1 = 0;
            f2 = 0;
            f3 = 0;
        end
        

        tspan = [time time+STEP];

        [~, qtemp] = ode45(@(time,y) qdotfunc(y, u, r), tspan, q);
        
        
        f = @(q, u) A(q(1),q(2),q(3),I1,I2,I3, tau1, tau2, tau3, f1, f2, f3, M,l,r) \ ... 
            b(q(1),q(2),q(3),u(1),u(2),u(3),I1,I2,I3, tau1, tau2, tau3, f1, f2, f3, M,l,r,g);
        
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

        A_bz = H * [0; 0; r/2];
        A_bx = H * [r/2; 0; 0];
        A_by = H * [0; r/2; 0];


       % bz_unit = bz ./ norm(bz);
        
        center(1) = q(4);
        center(2) = q(5);
        center(3) = r;
        %cg3arr = [cg3arr cg(3)];
        %p_ag = [cg(1) - q(4); cg(2) - q(5); cg(3)];

       % disp(q(3))
        disp(strcat('Time: ', num2str(time)));
        disp(strcat('X Position: ', num2str(q(4))));
        disp(strcat('Y Position: ', num2str(q(5))))


        plot3([center(1) center(1)+A_bz(1)], [center(2) center(2)+A_bz(2)], [center(3) center(3)+A_bz(3)], 'k');
        hold on;
        xlim([-3 3]);
        ylim([-3 3]);
        zlim([0 0.5]);
        grid ON;
        plot3([center(1) center(1)+A_bx(1)], [center(2) center(2)+A_bx(2)], [center(3) center(3)+A_bx(3)], 'r');
        plot3([center(1) center(1)+A_by(1)], [center(2) center(2)+A_by(2)], [center(3) center(3)+A_by(3)], 'g');
        %plot3([cg(1) cg(1)+bx(1)], [cg(2) cg(2)+bx(2)], [cg(3) cg(3)+bx(3)], 'r');
        plot3([q(4) center(1)],[q(5) center(2)],[0 center(3)], '-m');
        hold off;
        pause(0.0005);

        time = time + STEP;
        ind = ind + 1;

        %disp(time);

    end
end