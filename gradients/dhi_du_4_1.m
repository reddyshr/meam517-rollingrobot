function out1 = dhi_du_4_1(q1,q2,q3,q4,q5,u1,u2,u3,tau1,tau2,tau3,f1,f2,f3,q1p1,q2p1,q3p1,q4p1,q5p1,u1p1,u2p1,u3p1,tau1p1,tau2p1,tau3p1,f1p1,f2p1,f3p1,dt,M,l,r,g,I1,I2,I3)
%DHI_DU_4_1
%    OUT1 = DHI_DU_4_1(Q1,Q2,Q3,Q4,Q5,U1,U2,U3,TAU1,TAU2,TAU3,F1,F2,F3,Q1P1,Q2P1,Q3P1,Q4P1,Q5P1,U1P1,U2P1,U3P1,TAU1P1,TAU2P1,TAU3P1,F1P1,F2P1,F3P1,DT,M,L,R,G,I1,I2,I3)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    23-Nov-2018 20:52:54

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
t37 = t9.*t13;
t19 = t17-t37;
t20 = dt.*t19.*(1.0./8.0);
t21 = t14+t15+t20;
t22 = q3.*(1.0./2.0);
t23 = q3p1.*(1.0./2.0);
t24 = sin(q2);
t25 = sin(q2p1);
t26 = t9.*t13.*t25;
t39 = t3.*t7.*t24;
t27 = t26-t39+u3-u3p1;
t28 = dt.*t27.*(1.0./8.0);
t29 = t22+t23+t28;
t30 = l.^2;
t31 = M.^2;
t32 = r.^2;
t33 = t2.^2;
t34 = t4.^2;
t35 = t30.^2;
t36 = t32.^2;
t38 = sin(t21);
t40 = cos(t29);
t41 = cos(t21);
t42 = sin(t29);
t43 = q2.*(1.0./2.0);
t44 = q2p1.*(1.0./2.0);
t45 = t4.*u2;
t46 = t6.*u1;
t75 = t10.*u2p1;
t76 = t12.*u1p1;
t47 = t45+t46-t75-t76;
t48 = dt.*t47.*(1.0./8.0);
t49 = t43+t44+t48;
t50 = sin(t49);
t51 = I1.*t31.*t35;
t52 = I3.*t31.*t36;
t53 = I1.*I2.*I3;
t54 = I1.*I2.*M.*t30;
t55 = I1.*I3.*M.*t30;
t56 = I1.*I3.*M.*t32;
t57 = I2.*I3.*M.*t32;
t58 = I2.*t31.*t33.*t36;
t59 = I1.*t30.*t31.*t32;
t60 = I3.*t30.*t31.*t32;
t61 = I2.*t30.*t31.*t32.*t33;
t62 = I1.*I2.*M.*t32.*t33;
t63 = I1.*t31.*t33.*t34.*t36;
t64 = I1.*l.*r.*t2.*t4.*t31.*t32.*2.0;
t65 = I1.*l.*r.*t2.*t4.*t30.*t31.*4.0;
t66 = I3.*l.*r.*t2.*t4.*t31.*t32.*2.0;
t67 = I1.*t30.*t31.*t32.*t33.*t34.*5.0;
t68 = I1.*I3.*M.*t32.*t33.*t34;
t69 = I2.*l.*r.*t2.*t4.*t31.*t32.*t33.*2.0;
t70 = I1.*l.*r.*t2.*t4.*t31.*t32.*t33.*t34.*2.0;
t71 = I1.*I2.*M.*l.*r.*t2.*t4.*2.0;
t72 = I1.*I3.*M.*l.*r.*t2.*t4.*2.0;
t77 = I3.*t31.*t33.*t36;
t78 = I3.*t30.*t31.*t32.*t33;
t79 = I1.*I3.*M.*t32.*t33;
t80 = I2.*t31.*t33.*t34.*t36;
t81 = I2.*t30.*t31.*t32.*t33.*t34;
t82 = I2.*I3.*M.*t32.*t33.*t34;
t83 = I3.*l.*r.*t2.*t4.*t31.*t32.*t33.*2.0;
t84 = I2.*l.*r.*t2.*t4.*t31.*t32.*t33.*t34.*2.0;
t73 = t51+t52+t53+t54+t55+t56+t57+t58+t59+t60+t61+t62+t63+t64+t65+t66+t67+t68+t69+t70+t71+t72-t77-t78-t79-t80-t81-t82-t83-t84;
t74 = 1.0./t73;
out1 = -r.*(dt.*t74.*(t38.*t40+t41.*t42.*t50).*(t4.*t6.*t31.*t33.*t36+l.*r.*t2.*t6.*t30.*t31+l.*r.*t2.*t6.*t31.*t32+t4.*t6.*t30.*t31.*t32.*t33.*3.0+I3.*M.*l.*r.*t2.*t6+I3.*M.*t4.*t6.*t32.*t33+l.*r.*t2.*t6.*t31.*t32.*t33.*t34.*2.0).*(-1.0./8.0)+dt.*t74.*(t38.*t42-t40.*t41.*t50).*(t31.*t35+I2.*I3+I2.*M.*t30+I3.*M.*t30+I3.*M.*t32+t30.*t31.*t32+t31.*t33.*t34.*t36+I2.*M.*t32.*t33-I3.*M.*t32.*t33+I3.*M.*t32.*t33.*t34+t30.*t31.*t32.*t33.*t34.*5.0+l.*r.*t2.*t4.*t30.*t31.*4.0+l.*r.*t2.*t4.*t31.*t32.*2.0+I2.*M.*l.*r.*t2.*t4.*2.0+I3.*M.*l.*r.*t2.*t4.*2.0+l.*r.*t2.*t4.*t31.*t32.*t33.*t34.*2.0).*(1.0./8.0)+dt.*t41.*t74.*cos(t49).*(I2.*M.*l.*r.*t24+l.*r.*t24.*t30.*t31+l.*r.*t24.*t31.*t32+t2.*t4.*t24.*t31.*t36+t2.*t4.*t24.*t30.*t31.*t32.*3.0+I2.*M.*t2.*t4.*t24.*t32+l.*r.*t24.*t31.*t32.*t33.*t34.*2.0).*(1.0./8.0));
