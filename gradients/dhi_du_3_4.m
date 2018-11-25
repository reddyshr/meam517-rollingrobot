function out1 = dhi_du_3_4(q1,q2,q3,q4,q5,u1,u2,u3,tau1,tau2,tau3,f1,f2,f3,q1p1,q2p1,q3p1,q4p1,q5p1,u1p1,u2p1,u3p1,tau1p1,tau2p1,tau3p1,f1p1,f2p1,f3p1,dt,M,l,r,g,I1,I2,I3)
%DHI_DU_3_4
%    OUT1 = DHI_DU_3_4(Q1,Q2,Q3,Q4,Q5,U1,U2,U3,TAU1,TAU2,TAU3,F1,F2,F3,Q1P1,Q2P1,Q3P1,Q4P1,Q5P1,U1P1,U2P1,U3P1,TAU1P1,TAU2P1,TAU3P1,F1P1,F2P1,F3P1,DT,M,L,R,G,I1,I2,I3)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    23-Nov-2018 20:47:37

t2 = q2.*(1.0./2.0);
t3 = q2p1.*(1.0./2.0);
t4 = cos(q3);
t5 = t4.*u2;
t6 = cos(q3p1);
t7 = sin(q3);
t8 = t7.*u1;
t9 = sin(q3p1);
t10 = t5+t8-t9.*u1p1-t6.*u2p1;
t11 = dt.*t10.*(1.0./8.0);
t12 = t2+t3+t11;
t13 = r.^2;
t14 = cos(q2);
t15 = sin(q2);
t16 = t14.^2;
t17 = l.^2;
t18 = M.^2;
t19 = t13.^2;
t20 = t4.^2;
t21 = q3.*(1.0./2.0);
t22 = q3p1.*(1.0./2.0);
t23 = 1.0./t14;
t24 = t4.*u1;
t25 = t24-t7.*u2;
t26 = cos(q2p1);
t27 = 1.0./t26;
t28 = sin(q2p1);
t29 = t6.*u1p1;
t30 = t29-t9.*u2p1;
t31 = t27.*t28.*t30;
t32 = t31+u3-u3p1-t15.*t23.*t25;
t33 = dt.*t32.*(1.0./8.0);
t34 = t21+t22+t33;
t35 = t17.^2;
t36 = I1.*t18.*t35;
t37 = I3.*t18.*t19;
t38 = I1.*I2.*I3;
t39 = I1.*I2.*M.*t17;
t40 = I1.*I3.*M.*t17;
t41 = I1.*I3.*M.*t13;
t42 = I2.*I3.*M.*t13;
t43 = I2.*t16.*t18.*t19;
t44 = I1.*t13.*t17.*t18;
t45 = I3.*t13.*t17.*t18;
t46 = I2.*t13.*t16.*t17.*t18;
t47 = I1.*I2.*M.*t13.*t16;
t48 = I1.*t16.*t18.*t19.*t20;
t49 = I1.*l.*r.*t4.*t13.*t14.*t18.*2.0;
t50 = I1.*l.*r.*t4.*t14.*t17.*t18.*4.0;
t51 = I3.*l.*r.*t4.*t13.*t14.*t18.*2.0;
t52 = I1.*t13.*t16.*t17.*t18.*t20.*5.0;
t53 = I1.*I3.*M.*t13.*t16.*t20;
t54 = I2.*l.*r.*t4.*t13.*t14.*t16.*t18.*2.0;
t55 = I1.*l.*r.*t4.*t13.*t14.*t16.*t18.*t20.*2.0;
t56 = I1.*I2.*M.*l.*r.*t4.*t14.*2.0;
t57 = I1.*I3.*M.*l.*r.*t4.*t14.*2.0;
t60 = I3.*t16.*t18.*t19;
t61 = I3.*t13.*t16.*t17.*t18;
t62 = I1.*I3.*M.*t13.*t16;
t63 = I2.*t16.*t18.*t19.*t20;
t64 = I2.*t13.*t16.*t17.*t18.*t20;
t65 = I2.*I3.*M.*t13.*t16.*t20;
t66 = I3.*l.*r.*t4.*t13.*t14.*t16.*t18.*2.0;
t67 = I2.*l.*r.*t4.*t13.*t14.*t16.*t18.*t20.*2.0;
t58 = t36+t37+t38+t39+t40+t41+t42+t43+t44+t45+t46+t47+t48+t49+t50+t51+t52+t53+t54+t55+t56+t57-t60-t61-t62-t63-t64-t65-t66-t67;
t59 = 1.0./t58;
out1 = (sin(t12).*(dt.*t59.*cos(t34).*(I2.*M.*l.*t7.*t13.*t14.*t15-I3.*M.*l.*t7.*t13.*t14.*t15+I2.*M.*r.*t4.*t7.*t13.*t15.*t16-I3.*M.*r.*t4.*t7.*t13.*t15.*t16).*(1.0./8.0)-dt.*t59.*sin(t34).*(I1.*I3.*r.*t15+I3.*M.*r.*t13.*t15+I1.*M.*r.*t15.*t17+I1.*M.*l.*t4.*t13.*t14.*t15.*2.0+I1.*M.*r.*t13.*t15.*t16.*t20-I3.*M.*r.*t13.*t15.*t16.*t20).*(1.0./8.0)))./cos(t12)-dt.*t59.*(I1.*I2.*r.*t7.*t14+I2.*M.*r.*t7.*t13.*t14+I1.*M.*r.*t7.*t14.*t17+I1.*M.*l.*t4.*t7.*t13.*t16.*2.0+I1.*M.*r.*t7.*t13.*t14.*t16.*t20-I2.*M.*r.*t7.*t13.*t14.*t16.*t20).*(1.0./8.0);
