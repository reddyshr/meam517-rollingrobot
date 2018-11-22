function testsim(A,b)
        syms t q1 q2 q3 u1 u2 u3 u1_dot u2_dot u3_dot x_dot y_dot  l r M g tau f I1 I2 I3 real


    Mass = 1;
    I = 1.25*10^-3;
    J = 6.25*10^-4;
    radius = 0.02;
    h = 0.05;
    len = 0.03;
    grav = 9.81;

    theta = -pi/6;
    p = 72.84;
    s = 50;

    % Generalized Coordinates
    %q4 = x, q5 = y

    q = zeros(5,1);
    q(1) = theta/4;
    q(2) = theta/2;
    q(3) = theta;

    % Generalized Speeds

    u = zeros(5,1);
    u(1) = s+p*cos(theta);
    u(3) = 2*p*sin(theta);
    u(2) = p*sin(theta);

    s2 = sin(q(2));
    s3 = sin(q(3));

    c2 = cos(q(2));
    c3 = cos(q(3));
    
    u_1 = u(1);
    u_2 = u(2);
    u_3 = u(3);
    
    p = zeros(3,3);
    
    p(1,1) = (radius^2)*(s2^2 + (c2^2)*(s3^2)) + I/Mass;
    p(1,2) = radius*c2*s3*(h-radius + radius*c2*c3);
    p(1,3) = -radius*s2*(h-radius + radius*c2*c3);
    p(2,1) = p(1,2);
    p(2,2) = (radius^2)*(s2^2)+(h-radius+radius*c2*c3)^2 + J/Mass;
    p(2,3) = (radius^2)*s2*c2*s3;
    p(3,1) = p(1,3);
    p(3,2) = p(2,3);
    p(3,3) = (radius^2)*(c2^2)*(s3^2)+(h-radius+radius*c2*c3)^2 + J/Mass;
    
    n = zeros(3,1);
    n(1) = radius*(h-radius)*u_1*(s2*u_2+c2*s3*u_3);
    n(2) = (h-radius)*(radius*s2*((u_2^2) + (u_3^2)) + ...
            + (h-radius+radius*c2*c3)*u_3*u_1 + grav*s2) + u_3*u_1*(J-I)/Mass;
    n(3) = (h-radius)*(radius*c2*s3*((u_2^2) + (u_3^2)) + ...
        - (h-radius+radius*c2*c3)*u_1*u_2 + grav*c2*s3) + u_1*u_2*(I-J)/Mass;
    
    udot = p \ n;
    
    %My shit
    
    A = double(subs(A, [q1 q2 q3 u1 u2 u3 r l M g I1 I2 I3], [q(1) q(2) q(3) u_1 u_2 u_3 radius len Mass grav I J J]));
    b = double(subs(b, [q1 q2 q3 u1 u2 u3 r l M g I1 I2 I3], [q(1) q(2) q(3) u_1 u_2 u_3 radius len Mass grav I J J]));
    
    udot2 = A\b;

end