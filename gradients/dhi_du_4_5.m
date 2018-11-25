function out1 = dhi_du_4_5(q1,q2,q3,q4,q5,u1,u2,u3,tau1,tau2,tau3,f1,f2,f3,q1p1,q2p1,q3p1,q4p1,q5p1,u1p1,u2p1,u3p1,tau1p1,tau2p1,tau3p1,f1p1,f2p1,f3p1,dt,M,l,r,g,I1,I2,I3)
%DHI_DU_4_5
%    OUT1 = DHI_DU_4_5(Q1,Q2,Q3,Q4,Q5,U1,U2,U3,TAU1,TAU2,TAU3,F1,F2,F3,Q1P1,Q2P1,Q3P1,Q4P1,Q5P1,U1P1,U2P1,U3P1,TAU1P1,TAU2P1,TAU3P1,F1P1,F2P1,F3P1,DT,M,L,R,G,I1,I2,I3)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    23-Nov-2018 20:57:52

t2 = cos(q2);
t3 = 1.0./t2;
t4 = cos(q3);
t5 = t4.*u1;
t6 = sin(q3);
t16 = t6.*u2;
t7 = t5-t16;
t8 = cos(q2p1);
t9 = 1.0./t8;
t10 = cos(q3p1);
t11 = t10.*u1p1;
t12 = sin(q3p1);
t18 = t12.*u2p1;
t13 = t11-t18;
t14 = q1.*(1.0./2.0);
t15 = q1p1.*(1.0./2.0);
t17 = t3.*t7;
t36 = t9.*t13;
t19 = t17-t36;
t20 = dt.*t19.*(1.0./8.0);
t21 = t14+t15+t20;
t22 = q3.*(1.0./2.0);
t23 = q3p1.*(1.0./2.0);
t24 = sin(q2);
t25 = sin(q2p1);
t26 = t9.*t13.*t25;
t38 = t3.*t7.*t24;
t27 = t26-t38+u3-u3p1;
t28 = dt.*t27.*(1.0./8.0);
t29 = t22+t23+t28;
t30 = r.^2;
t31 = t2.^2;
t32 = t4.^2;
t33 = l.^2;
t34 = M.^2;
t35 = t30.^2;
t37 = sin(t21);
t39 = cos(t29);
t40 = cos(t21);
t41 = sin(t29);
t42 = q2.*(1.0./2.0);
t43 = q2p1.*(1.0./2.0);
t44 = t4.*u2;
t45 = t6.*u1;
t75 = t10.*u2p1;
t76 = t12.*u1p1;
t46 = t44+t45-t75-t76;
t47 = dt.*t46.*(1.0./8.0);
t48 = t42+t43+t47;
t49 = sin(t48);
t50 = t33.^2;
t51 = I1.*t34.*t50;
t52 = I3.*t34.*t35;
t53 = I1.*I2.*I3;
t54 = I1.*I2.*M.*t33;
t55 = I1.*I3.*M.*t33;
t56 = I1.*I3.*M.*t30;
t57 = I2.*I3.*M.*t30;
t58 = I2.*t31.*t34.*t35;
t59 = I1.*t30.*t33.*t34;
t60 = I3.*t30.*t33.*t34;
t61 = I2.*t30.*t31.*t33.*t34;
t62 = I1.*I2.*M.*t30.*t31;
t63 = I1.*t31.*t32.*t34.*t35;
t64 = I1.*l.*r.*t2.*t4.*t30.*t34.*2.0;
t65 = I1.*l.*r.*t2.*t4.*t33.*t34.*4.0;
t66 = I3.*l.*r.*t2.*t4.*t30.*t34.*2.0;
t67 = I1.*t30.*t31.*t32.*t33.*t34.*5.0;
t68 = I1.*I3.*M.*t30.*t31.*t32;
t69 = I2.*l.*r.*t2.*t4.*t30.*t31.*t34.*2.0;
t70 = I1.*l.*r.*t2.*t4.*t30.*t31.*t32.*t34.*2.0;
t71 = I1.*I2.*M.*l.*r.*t2.*t4.*2.0;
t72 = I1.*I3.*M.*l.*r.*t2.*t4.*2.0;
t77 = I3.*t31.*t34.*t35;
t78 = I3.*t30.*t31.*t33.*t34;
t79 = I1.*I3.*M.*t30.*t31;
t80 = I2.*t31.*t32.*t34.*t35;
t81 = I2.*t30.*t31.*t32.*t33.*t34;
t82 = I2.*I3.*M.*t30.*t31.*t32;
t83 = I3.*l.*r.*t2.*t4.*t30.*t31.*t34.*2.0;
t84 = I2.*l.*r.*t2.*t4.*t30.*t31.*t32.*t34.*2.0;
t73 = t51+t52+t53+t54+t55+t56+t57+t58+t59+t60+t61+t62+t63+t64+t65+t66+t67+t68+t69+t70+t71+t72-t77-t78-t79-t80-t81-t82-t83-t84;
t74 = 1.0./t73;
out1 = r.*(dt.*t74.*(t37.*t41-t39.*t40.*t49).*(I2.*I3.*r.*t24+I3.*M.*r.*t24.*t30+I3.*M.*r.*t24.*t33+I2.*M.*r.*t24.*t30.*t31-I3.*M.*r.*t24.*t30.*t31+I3.*M.*l.*t2.*t4.*t24.*t30.*2.0-I2.*M.*r.*t24.*t30.*t31.*t32+I3.*M.*r.*t24.*t30.*t31.*t32).*(1.0./8.0)+dt.*t74.*(t37.*t39+t40.*t41.*t49).*(I1.*M.*l.*t2.*t6.*t24.*t30-I3.*M.*l.*t2.*t6.*t24.*t30+I1.*M.*r.*t4.*t6.*t24.*t30.*t31-I3.*M.*r.*t4.*t6.*t24.*t30.*t31).*(1.0./8.0)-dt.*t40.*t74.*cos(t48).*(I1.*I2.*l+I1.*M.*l.*t30+I1.*M.*l.*t33-I1.*M.*l.*t30.*t31+I2.*M.*l.*t30.*t31+I1.*I2.*r.*t2.*t4+I1.*M.*l.*t30.*t31.*t32.*3.0-I2.*M.*l.*t30.*t31.*t32+I1.*M.*r.*t2.*t4.*t30+I1.*M.*r.*t2.*t4.*t33.*3.0-I1.*M.*r.*t2.*t4.*t30.*t31+I2.*M.*r.*t2.*t4.*t30.*t31+I1.*M.*r.*t2.*t4.*t30.*t31.*t32-I2.*M.*r.*t2.*t4.*t30.*t31.*t32).*(1.0./8.0));
