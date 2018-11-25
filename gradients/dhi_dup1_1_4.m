function out1 = dhi_dup1_1_4(q1,q2,q3,q4,q5,u1,u2,u3,tau1,tau2,tau3,f1,f2,f3,q1p1,q2p1,q3p1,q4p1,q5p1,u1p1,u2p1,u3p1,tau1p1,tau2p1,tau3p1,f1p1,f2p1,f3p1,dt,M,l,r,g,I1,I2,I3)
%DHI_DUP1_1_4
%    OUT1 = DHI_DUP1_1_4(Q1,Q2,Q3,Q4,Q5,U1,U2,U3,TAU1,TAU2,TAU3,F1,F2,F3,Q1P1,Q2P1,Q3P1,Q4P1,Q5P1,U1P1,U2P1,U3P1,TAU1P1,TAU2P1,TAU3P1,F1P1,F2P1,F3P1,DT,M,L,R,G,I1,I2,I3)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    23-Nov-2018 20:35:37

t2 = cos(q3);
t3 = sin(q3);
t4 = cos(q3p1);
t5 = sin(q3p1);
t6 = r.^2;
t7 = cos(q2p1);
t8 = sin(q2p1);
t9 = t7.^2;
t10 = l.^2;
t11 = M.^2;
t12 = t6.^2;
t13 = t4.^2;
t14 = q3.*(1.0./2.0);
t15 = q3p1.*(1.0./2.0);
t16 = cos(q2);
t17 = 1.0./t16;
t18 = sin(q2);
t19 = t2.*u1;
t20 = t19-t3.*u2;
t21 = 1.0./t7;
t22 = t4.*u1p1;
t23 = t22-t5.*u2p1;
t24 = t8.*t21.*t23;
t25 = t24+u3-u3p1-t17.*t18.*t20;
t26 = dt.*t25.*(1.0./8.0);
t27 = t14+t15+t26;
t28 = t10.^2;
t29 = I1.*t11.*t28;
t30 = I3.*t11.*t12;
t31 = I1.*I2.*I3;
t32 = I1.*I2.*M.*t10;
t33 = I1.*I3.*M.*t10;
t34 = I1.*I3.*M.*t6;
t35 = I2.*I3.*M.*t6;
t36 = I2.*t9.*t11.*t12;
t37 = I1.*t6.*t10.*t11;
t38 = I3.*t6.*t10.*t11;
t39 = I2.*t6.*t9.*t10.*t11;
t40 = I1.*I2.*M.*t6.*t9;
t41 = I1.*t9.*t11.*t12.*t13;
t42 = I1.*l.*r.*t4.*t6.*t7.*t11.*2.0;
t43 = I1.*l.*r.*t4.*t7.*t10.*t11.*4.0;
t44 = I3.*l.*r.*t4.*t6.*t7.*t11.*2.0;
t45 = I1.*t6.*t9.*t10.*t11.*t13.*5.0;
t46 = I1.*I3.*M.*t6.*t9.*t13;
t47 = I2.*l.*r.*t4.*t6.*t7.*t9.*t11.*2.0;
t48 = I1.*l.*r.*t4.*t6.*t7.*t9.*t11.*t13.*2.0;
t49 = I1.*I2.*M.*l.*r.*t4.*t7.*2.0;
t50 = I1.*I3.*M.*l.*r.*t4.*t7.*2.0;
t51 = t29+t30+t31+t32+t33+t34+t35+t36+t37+t38+t39+t40+t41+t42+t43+t44+t45+t46+t47+t48+t49+t50-I3.*t9.*t11.*t12-I1.*I3.*M.*t6.*t9-I3.*t6.*t9.*t10.*t11-I2.*t9.*t11.*t12.*t13-I2.*I3.*M.*t6.*t9.*t13-I2.*t6.*t9.*t10.*t11.*t13-I3.*l.*r.*t4.*t6.*t7.*t9.*t11.*2.0-I2.*l.*r.*t4.*t6.*t7.*t9.*t11.*t13.*2.0;
t52 = 1.0./t51;
out1 = (dt.*t52.*cos(t27).*(I2.*M.*l.*t5.*t6.*t7.*t8-I3.*M.*l.*t5.*t6.*t7.*t8+I2.*M.*r.*t4.*t5.*t6.*t8.*t9-I3.*M.*r.*t4.*t5.*t6.*t8.*t9).*(1.0./8.0)-dt.*t52.*sin(t27).*(I1.*I3.*r.*t8+I3.*M.*r.*t6.*t8+I1.*M.*r.*t8.*t10+I1.*M.*l.*t4.*t6.*t7.*t8.*2.0+I1.*M.*r.*t6.*t8.*t9.*t13-I3.*M.*r.*t6.*t8.*t9.*t13).*(1.0./8.0))./cos(q2.*(1.0./2.0)+q2p1.*(1.0./2.0)+dt.*(t2.*u2+t3.*u1-t5.*u1p1-t4.*u2p1).*(1.0./8.0));
