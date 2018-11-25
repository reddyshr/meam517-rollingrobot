function out1 = dhi_dxp1_4_1(q1,q2,q3,q4,q5,u1,u2,u3,tau1,tau2,tau3,f1,f2,f3,q1p1,q2p1,q3p1,q4p1,q5p1,u1p1,u2p1,u3p1,tau1p1,tau2p1,tau3p1,f1p1,f2p1,f3p1,dt,M,l,r,g,I1,I2,I3)
%DHI_DXP1_4_1
%    OUT1 = DHI_DXP1_4_1(Q1,Q2,Q3,Q4,Q5,U1,U2,U3,TAU1,TAU2,TAU3,F1,F2,F3,Q1P1,Q2P1,Q3P1,Q4P1,Q5P1,U1P1,U2P1,U3P1,TAU1P1,TAU2P1,TAU3P1,F1P1,F2P1,F3P1,DT,M,L,R,G,I1,I2,I3)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    23-Nov-2018 20:53:19

t2 = cos(q1p1);
t3 = cos(q3p1);
t4 = sin(q1p1);
t5 = sin(q2p1);
t6 = sin(q3p1);
t7 = cos(q2p1);
t8 = cos(q2);
t9 = 1.0./t8;
t10 = cos(q3);
t11 = t10.*u1;
t12 = sin(q3);
t19 = t12.*u2;
t13 = t11-t19;
t14 = 1.0./t7;
t15 = t3.*u1p1;
t21 = t6.*u2p1;
t16 = t15-t21;
t17 = q1.*(1.0./2.0);
t18 = q1p1.*(1.0./2.0);
t20 = t9.*t13;
t56 = t14.*t16;
t22 = t20-t56;
t23 = dt.*t22.*(1.0./8.0);
t24 = t17+t18+t23;
t25 = q3.*(1.0./2.0);
t26 = q3p1.*(1.0./2.0);
t27 = sin(q2);
t28 = t5.*t14.*t16;
t58 = t9.*t13.*t27;
t29 = t28-t58+u3-u3p1;
t30 = dt.*t29.*(1.0./8.0);
t31 = t25+t26+t30;
t32 = l.^2;
t33 = M.^2;
t34 = r.^2;
t35 = t34.^2;
t36 = t8.^2;
t37 = t10.^2;
t38 = I1.^2;
t39 = t32.^2;
t40 = I3.^2;
t41 = u2.^2;
t42 = u3.^2;
t43 = I1.*t33.*t39;
t44 = I3.*t33.*t35;
t45 = I1.*I2.*I3;
t46 = I1.*I2.*M.*t32;
t47 = I1.*I3.*M.*t32;
t48 = I1.*I3.*M.*t34;
t49 = I2.*I3.*M.*t34;
t50 = t7.^2;
t51 = I1.*t32.*t33.*t34;
t52 = I3.*t32.*t33.*t34;
t53 = t3.^2;
t54 = u2p1.^2;
t55 = u3p1.^2;
t57 = cos(t24);
t59 = sin(t31);
t60 = sin(t24);
t61 = cos(t31);
t62 = q2.*(1.0./2.0);
t63 = q2p1.*(1.0./2.0);
t64 = t10.*u2;
t65 = t12.*u1;
t101 = t3.*u2p1;
t102 = t6.*u1p1;
t66 = t64+t65-t101-t102;
t67 = dt.*t66.*(1.0./8.0);
t68 = t62+t63+t67;
t69 = sin(t68);
t70 = I2.*t33.*t35.*t36;
t71 = I2.*t32.*t33.*t34.*t36;
t72 = I1.*I2.*M.*t34.*t36;
t73 = I1.*t33.*t35.*t36.*t37;
t74 = I1.*l.*r.*t8.*t10.*t33.*t34.*2.0;
t75 = I1.*l.*r.*t8.*t10.*t32.*t33.*4.0;
t76 = I3.*l.*r.*t8.*t10.*t33.*t34.*2.0;
t77 = I1.*t32.*t33.*t34.*t36.*t37.*5.0;
t78 = I1.*I3.*M.*t34.*t36.*t37;
t79 = I2.*l.*r.*t8.*t10.*t33.*t34.*t36.*2.0;
t80 = I1.*l.*r.*t8.*t10.*t33.*t34.*t36.*t37.*2.0;
t81 = I1.*I2.*M.*l.*r.*t8.*t10.*2.0;
t82 = I1.*I3.*M.*l.*r.*t8.*t10.*2.0;
t103 = I3.*t33.*t35.*t36;
t104 = I3.*t32.*t33.*t34.*t36;
t105 = I1.*I3.*M.*t34.*t36;
t106 = I2.*t33.*t35.*t36.*t37;
t107 = I2.*t32.*t33.*t34.*t36.*t37;
t108 = I2.*I3.*M.*t34.*t36.*t37;
t109 = I3.*l.*r.*t8.*t10.*t33.*t34.*t36.*2.0;
t110 = I2.*l.*r.*t8.*t10.*t33.*t34.*t36.*t37.*2.0;
t83 = t43+t44+t45+t46+t47+t48+t49+t51+t52+t70+t71+t72+t73+t74+t75+t76+t77+t78+t79+t80+t81+t82-t103-t104-t105-t106-t107-t108-t109-t110;
t84 = 1.0./t83;
t85 = I2.^2;
t86 = I2.*t33.*t35.*t50;
t87 = I2.*t32.*t33.*t34.*t50;
t88 = I1.*I2.*M.*t34.*t50;
t89 = I1.*t33.*t35.*t50.*t53;
t90 = I1.*l.*r.*t3.*t7.*t33.*t34.*2.0;
t91 = I1.*l.*r.*t3.*t7.*t32.*t33.*4.0;
t92 = I3.*l.*r.*t3.*t7.*t33.*t34.*2.0;
t93 = I1.*t32.*t33.*t34.*t50.*t53.*5.0;
t94 = I1.*I3.*M.*t34.*t50.*t53;
t95 = I2.*l.*r.*t3.*t7.*t33.*t34.*t50.*2.0;
t96 = I1.*l.*r.*t3.*t7.*t33.*t34.*t50.*t53.*2.0;
t97 = I1.*I2.*M.*l.*r.*t3.*t7.*2.0;
t98 = I1.*I3.*M.*l.*r.*t3.*t7.*2.0;
t111 = I3.*t33.*t35.*t50;
t112 = I3.*t32.*t33.*t34.*t50;
t113 = I1.*I3.*M.*t34.*t50;
t114 = I2.*t33.*t35.*t50.*t53;
t115 = I2.*t32.*t33.*t34.*t50.*t53;
t116 = I2.*I3.*M.*t34.*t50.*t53;
t117 = I3.*l.*r.*t3.*t7.*t33.*t34.*t50.*2.0;
t118 = I2.*l.*r.*t3.*t7.*t33.*t34.*t50.*t53.*2.0;
t99 = t43+t44+t45+t46+t47+t48+t49+t51+t52+t86+t87+t88+t89+t90+t91+t92+t93+t94+t95+t96+t97+t98-t111-t112-t113-t114-t115-t116-t117-t118;
t100 = 1.0./t99;
out1 = -r.*((t57.*t61.*(1.0./2.0)-t59.*t60.*t69.*(1.0./2.0)).*(u2.*(1.0./2.0)+u2p1.*(1.0./2.0)+dt.*(t84.*(I1.*I3.*tau2+I1.*t40.*u1.*u3-I3.*t38.*u1.*u3+t33.*t35.*t36.*tau2-I1.*I3.*f3.*l+I1.*M.*t32.*tau2+I3.*M.*t34.*tau2-I1.*M.*f3.*l.*t32-I3.*M.*f3.*l.*t34+I1.*I3.*f1.*r.*t27+I1.*M.*t34.*t36.*tau2+I1.*t33.*t39.*u1.*u3-M.*t32.*t38.*u1.*u3+M.*t34.*t40.*u1.*u3+t32.*t33.*t34.*t36.*tau2-t33.*t35.*t36.*t37.*tau2-t8.*t12.*t27.*t33.*t35.*tau3-t10.*t12.*t33.*t35.*t36.*tau1-t32.*t33.*t34.*t36.*t37.*tau2+I1.*I3.*M.*g.*l.*t27+I1.*I3.*M.*t32.*u1.*u3.*2.0-I1.*I3.*M.*t34.*u1.*u3-I1.*M.*f3.*l.*t34.*t36+I3.*M.*f3.*l.*t34.*t36-I1.*I3.*f3.*r.*t8.*t10+I1.*M.*f1.*r.*t27.*t32+I3.*M.*f1.*r.*t27.*t34-I3.*M.*t34.*t36.*t37.*tau2+I1.*g.*l.*t27.*t32.*t33+I3.*g.*l.*t27.*t33.*t34+I3.*t32.*t33.*t34.*u1.*u3-I1.*t33.*t35.*t36.*u1.*u3+I3.*t33.*t35.*t36.*u1.*u3-M.*t34.*t36.*t38.*u1.*u3+I1.*I3.*M.*l.*r.*t27.*t41+I1.*I3.*M.*l.*r.*t27.*t42+I1.*I3.*M.*t34.*t36.*u1.*u3-I1.*M.*f3.*l.*t34.*t36.*t37.*2.0-I1.*M.*f3.*r.*t8.*t10.*t32.*3.0-I3.*M.*f3.*r.*t8.*t10.*t34+I1.*M.*l.*r.*t8.*t10.*tau2.*2.0-I3.*M.*l.*r.*t8.*t12.*tau1-I1.*M.*t8.*t12.*t27.*t34.*tau3-I3.*M.*t10.*t12.*t34.*t36.*tau1+I1.*l.*r.*t27.*t32.*t33.*t41+I1.*l.*r.*t27.*t32.*t33.*t42+I3.*l.*r.*t27.*t33.*t34.*t41+I3.*l.*r.*t27.*t33.*t34.*t42+I1.*t33.*t35.*t36.*t37.*u1.*u3-I3.*t33.*t35.*t36.*t37.*u1.*u3-M.*t34.*t36.*t37.*t40.*u1.*u3-l.*r.*t8.*t12.*t32.*t33.*tau1-l.*r.*t8.*t12.*t33.*t34.*tau1-t8.*t12.*t27.*t32.*t33.*t34.*tau3-t10.*t12.*t32.*t33.*t34.*t36.*tau1.*3.0+I1.*I3.*M.*t34.*t36.*t37.*u1.*u3+I1.*M.*f1.*l.*t8.*t10.*t27.*t34.*2.0-I1.*M.*f2.*l.*t8.*t12.*t27.*t34+I3.*M.*f2.*l.*t8.*t12.*t27.*t34-I1.*M.*f3.*r.*t8.*t10.*t34.*t36+I3.*M.*f3.*r.*t8.*t10.*t34.*t36+I1.*M.*f1.*r.*t27.*t34.*t36.*t37-I3.*M.*f1.*r.*t27.*t34.*t36.*t37+I1.*g.*l.*t27.*t33.*t34.*t36.*t37-I3.*g.*l.*t27.*t33.*t34.*t36.*t37+I1.*g.*r.*t8.*t10.*t27.*t32.*t33.*2.0-M.*l.*r.*t8.*t10.*t38.*u1.*u3.*2.0+M.*l.*r.*t8.*t12.*t40.*u2.*u3+I1.*t8.*t10.*t27.*t32.*t33.*t34.*t41.*2.0+I1.*t8.*t10.*t27.*t32.*t33.*t34.*t42.*2.0-I1.*t8.*t12.*t27.*t33.*t35.*u1.*u2+I2.*t8.*t12.*t27.*t33.*t35.*u1.*u2-I2.*t10.*t12.*t33.*t35.*t36.*u2.*u3+I3.*t10.*t12.*t33.*t35.*t36.*u2.*u3+I1.*t32.*t33.*t34.*t36.*t37.*u1.*u3.*3.0-I3.*t32.*t33.*t34.*t36.*t37.*u1.*u3-M.*t8.*t12.*t27.*t34.*t38.*u1.*u2+M.*t10.*t12.*t34.*t36.*t40.*u2.*u3+l.*r.*t8.*t10.*t33.*t34.*t36.*tau2.*2.0-I1.*M.*f2.*r.*t10.*t12.*t27.*t34.*t36+I3.*M.*f2.*r.*t10.*t12.*t27.*t34.*t36+I1.*l.*r.*t27.*t33.*t34.*t36.*t37.*t41+I1.*l.*r.*t27.*t33.*t34.*t36.*t37.*t42-I3.*l.*r.*t27.*t33.*t34.*t36.*t37.*t41-I3.*l.*r.*t27.*t33.*t34.*t36.*t37.*t42+I1.*l.*r.*t8.*t10.*t32.*t33.*u1.*u3.*3.0-I2.*l.*r.*t8.*t12.*t32.*t33.*u2.*u3+I3.*l.*r.*t8.*t10.*t33.*t34.*u1.*u3+I3.*l.*r.*t8.*t12.*t32.*t33.*u2.*u3-I2.*l.*r.*t8.*t12.*t33.*t34.*u2.*u3+I3.*l.*r.*t8.*t12.*t33.*t34.*u2.*u3+I2.*t8.*t12.*t27.*t32.*t33.*t34.*u1.*u2-I3.*t8.*t12.*t27.*t32.*t33.*t34.*u1.*u2-I2.*t10.*t12.*t32.*t33.*t34.*t36.*u2.*u3.*3.0+I3.*t10.*t12.*t32.*t33.*t34.*t36.*u2.*u3.*3.0-l.*r.*t10.*t12.*t27.*t33.*t34.*t36.*tau3.*2.0-l.*r.*t8.*t10.*t33.*t34.*t36.*t37.*tau2.*2.0-l.*r.*t8.*t12.*t33.*t34.*t36.*t37.*tau1.*2.0+I1.*I3.*M.*l.*r.*t8.*t10.*u1.*u3.*3.0-I2.*I3.*M.*l.*r.*t8.*t12.*u2.*u3+I1.*I2.*M.*t8.*t12.*t27.*t34.*u1.*u2-I2.*I3.*M.*t10.*t12.*t34.*t36.*u2.*u3-I1.*l.*r.*t8.*t10.*t33.*t34.*t36.*u1.*u3+I3.*l.*r.*t8.*t10.*t33.*t34.*t36.*u1.*u3-I1.*l.*r.*t10.*t12.*t27.*t33.*t34.*t36.*u1.*u2+I2.*l.*r.*t10.*t12.*t27.*t33.*t34.*t36.*u1.*u2.*2.0-I3.*l.*r.*t10.*t12.*t27.*t33.*t34.*t36.*u1.*u2+I1.*l.*r.*t8.*t10.*t33.*t34.*t36.*t37.*u1.*u3.*2.0-I3.*l.*r.*t8.*t10.*t33.*t34.*t36.*t37.*u1.*u3.*2.0-I2.*l.*r.*t8.*t12.*t33.*t34.*t36.*t37.*u2.*u3.*2.0+I3.*l.*r.*t8.*t12.*t33.*t34.*t36.*t37.*u2.*u3.*2.0)-t100.*(I1.*I3.*tau2p1+I1.*t40.*u1p1.*u3p1-I3.*t38.*u1p1.*u3p1+t33.*t35.*t50.*tau2p1-I1.*I3.*f3p1.*l+I1.*M.*t32.*tau2p1+I3.*M.*t34.*tau2p1-I1.*M.*f3p1.*l.*t32-I3.*M.*f3p1.*l.*t34+I1.*I3.*f1p1.*r.*t5+I1.*M.*t34.*t50.*tau2p1+I1.*t33.*t39.*u1p1.*u3p1-M.*t32.*t38.*u1p1.*u3p1+M.*t34.*t40.*u1p1.*u3p1+t32.*t33.*t34.*t50.*tau2p1-t33.*t35.*t50.*t53.*tau2p1-t5.*t6.*t7.*t33.*t35.*tau3p1-t3.*t6.*t33.*t35.*t50.*tau1p1-t32.*t33.*t34.*t50.*t53.*tau2p1+I1.*I3.*M.*g.*l.*t5+I1.*I3.*M.*t32.*u1p1.*u3p1.*2.0-I1.*I3.*M.*t34.*u1p1.*u3p1-I1.*M.*f3p1.*l.*t34.*t50+I3.*M.*f3p1.*l.*t34.*t50-I1.*I3.*f3p1.*r.*t3.*t7+I1.*M.*f1p1.*r.*t5.*t32+I3.*M.*f1p1.*r.*t5.*t34-I3.*M.*t34.*t50.*t53.*tau2p1+I1.*g.*l.*t5.*t32.*t33+I3.*g.*l.*t5.*t33.*t34+I3.*t32.*t33.*t34.*u1p1.*u3p1-I1.*t33.*t35.*t50.*u1p1.*u3p1+I3.*t33.*t35.*t50.*u1p1.*u3p1-M.*t34.*t38.*t50.*u1p1.*u3p1+I1.*I3.*M.*l.*r.*t5.*t54+I1.*I3.*M.*l.*r.*t5.*t55+I1.*I3.*M.*t34.*t50.*u1p1.*u3p1-I1.*M.*f3p1.*l.*t34.*t50.*t53.*2.0-I1.*M.*f3p1.*r.*t3.*t7.*t32.*3.0-I3.*M.*f3p1.*r.*t3.*t7.*t34-I3.*M.*l.*r.*t6.*t7.*tau1p1+I1.*M.*l.*r.*t3.*t7.*tau2p1.*2.0-I1.*M.*t5.*t6.*t7.*t34.*tau3p1-I3.*M.*t3.*t6.*t34.*t50.*tau1p1+I1.*l.*r.*t5.*t32.*t33.*t54+I1.*l.*r.*t5.*t32.*t33.*t55+I3.*l.*r.*t5.*t33.*t34.*t54+I3.*l.*r.*t5.*t33.*t34.*t55+I1.*t33.*t35.*t50.*t53.*u1p1.*u3p1-I3.*t33.*t35.*t50.*t53.*u1p1.*u3p1-M.*t34.*t40.*t50.*t53.*u1p1.*u3p1-l.*r.*t6.*t7.*t32.*t33.*tau1p1-l.*r.*t6.*t7.*t33.*t34.*tau1p1-t5.*t6.*t7.*t32.*t33.*t34.*tau3p1-t3.*t6.*t32.*t33.*t34.*t50.*tau1p1.*3.0+I1.*I3.*M.*t34.*t50.*t53.*u1p1.*u3p1+I1.*M.*f1p1.*l.*t3.*t5.*t7.*t34.*2.0-I1.*M.*f2p1.*l.*t5.*t6.*t7.*t34+I3.*M.*f2p1.*l.*t5.*t6.*t7.*t34-I1.*M.*f3p1.*r.*t3.*t7.*t34.*t50+I3.*M.*f3p1.*r.*t3.*t7.*t34.*t50+I1.*M.*f1p1.*r.*t5.*t34.*t50.*t53-I3.*M.*f1p1.*r.*t5.*t34.*t50.*t53+I1.*g.*l.*t5.*t33.*t34.*t50.*t53-I3.*g.*l.*t5.*t33.*t34.*t50.*t53+I1.*g.*r.*t3.*t5.*t7.*t32.*t33.*2.0-M.*l.*r.*t3.*t7.*t38.*u1p1.*u3p1.*2.0+M.*l.*r.*t6.*t7.*t40.*u2p1.*u3p1+I1.*t3.*t5.*t7.*t32.*t33.*t34.*t54.*2.0+I1.*t3.*t5.*t7.*t32.*t33.*t34.*t55.*2.0-I1.*t5.*t6.*t7.*t33.*t35.*u1p1.*u2p1+I2.*t5.*t6.*t7.*t33.*t35.*u1p1.*u2p1-I2.*t3.*t6.*t33.*t35.*t50.*u2p1.*u3p1+I3.*t3.*t6.*t33.*t35.*t50.*u2p1.*u3p1+I1.*t32.*t33.*t34.*t50.*t53.*u1p1.*u3p1.*3.0-I3.*t32.*t33.*t34.*t50.*t53.*u1p1.*u3p1-M.*t5.*t6.*t7.*t34.*t38.*u1p1.*u2p1+M.*t3.*t6.*t34.*t40.*t50.*u2p1.*u3p1+l.*r.*t3.*t7.*t33.*t34.*t50.*tau2p1.*2.0-I1.*M.*f2p1.*r.*t3.*t5.*t6.*t34.*t50+I3.*M.*f2p1.*r.*t3.*t5.*t6.*t34.*t50+I1.*l.*r.*t5.*t33.*t34.*t50.*t53.*t54+I1.*l.*r.*t5.*t33.*t34.*t50.*t53.*t55-I3.*l.*r.*t5.*t33.*t34.*t50.*t53.*t54-I3.*l.*r.*t5.*t33.*t34.*t50.*t53.*t55+I1.*l.*r.*t3.*t7.*t32.*t33.*u1p1.*u3p1.*3.0+I3.*l.*r.*t3.*t7.*t33.*t34.*u1p1.*u3p1-I2.*l.*r.*t6.*t7.*t32.*t33.*u2p1.*u3p1+I3.*l.*r.*t6.*t7.*t32.*t33.*u2p1.*u3p1-I2.*l.*r.*t6.*t7.*t33.*t34.*u2p1.*u3p1+I3.*l.*r.*t6.*t7.*t33.*t34.*u2p1.*u3p1+I2.*t5.*t6.*t7.*t32.*t33.*t34.*u1p1.*u2p1-I3.*t5.*t6.*t7.*t32.*t33.*t34.*u1p1.*u2p1-I2.*t3.*t6.*t32.*t33.*t34.*t50.*u2p1.*u3p1.*3.0+I3.*t3.*t6.*t32.*t33.*t34.*t50.*u2p1.*u3p1.*3.0-l.*r.*t3.*t5.*t6.*t33.*t34.*t50.*tau3p1.*2.0-l.*r.*t6.*t7.*t33.*t34.*t50.*t53.*tau1p1.*2.0-l.*r.*t3.*t7.*t33.*t34.*t50.*t53.*tau2p1.*2.0+I1.*I3.*M.*l.*r.*t3.*t7.*u1p1.*u3p1.*3.0-I2.*I3.*M.*l.*r.*t6.*t7.*u2p1.*u3p1+I1.*I2.*M.*t5.*t6.*t7.*t34.*u1p1.*u2p1-I2.*I3.*M.*t3.*t6.*t34.*t50.*u2p1.*u3p1-I1.*l.*r.*t3.*t7.*t33.*t34.*t50.*u1p1.*u3p1+I3.*l.*r.*t3.*t7.*t33.*t34.*t50.*u1p1.*u3p1-I1.*l.*r.*t3.*t5.*t6.*t33.*t34.*t50.*u1p1.*u2p1+I2.*l.*r.*t3.*t5.*t6.*t33.*t34.*t50.*u1p1.*u2p1.*2.0-I3.*l.*r.*t3.*t5.*t6.*t33.*t34.*t50.*u1p1.*u2p1+I1.*l.*r.*t3.*t7.*t33.*t34.*t50.*t53.*u1p1.*u3p1.*2.0-I3.*l.*r.*t3.*t7.*t33.*t34.*t50.*t53.*u1p1.*u3p1.*2.0-I2.*l.*r.*t6.*t7.*t33.*t34.*t50.*t53.*u2p1.*u3p1.*2.0+I3.*l.*r.*t6.*t7.*t33.*t34.*t50.*t53.*u2p1.*u3p1.*2.0)).*(1.0./8.0))+(t57.*t59.*(1.0./2.0)+t60.*t61.*t69.*(1.0./2.0)).*(u1.*(1.0./2.0)+u1p1.*(1.0./2.0)+dt.*(t84.*(I2.*I3.*tau1+t33.*t39.*tau1-I2.*t40.*u2.*u3+I3.*t85.*u2.*u3+t32.*t33.*t34.*tau1+I2.*M.*t32.*tau1+I3.*M.*t32.*tau1+I3.*M.*t34.*tau1-I2.*I3.*f2.*r.*t27+I2.*M.*t34.*t36.*tau1-I3.*M.*t34.*t36.*tau1+I2.*t33.*t39.*u2.*u3-I3.*t33.*t39.*u2.*u3-M.*t32.*t40.*u2.*u3-M.*t34.*t40.*u2.*u3+M.*t32.*t85.*u2.*u3+t33.*t35.*t36.*t37.*tau1+l.*r.*t27.*t32.*t33.*tau3+l.*r.*t27.*t33.*t34.*tau3+t8.*t10.*t27.*t33.*t35.*tau3-t10.*t12.*t33.*t35.*t36.*tau2+t32.*t33.*t34.*t36.*t37.*tau1.*5.0+I2.*I3.*M.*t34.*u2.*u3-I2.*I3.*f3.*r.*t8.*t12-I3.*M.*f2.*r.*t27.*t32-I3.*M.*f2.*r.*t27.*t34+I2.*M.*l.*r.*t27.*tau3+I3.*M.*t34.*t36.*t37.*tau1+I2.*t32.*t33.*t34.*u2.*u3-I3.*t32.*t33.*t34.*u2.*u3+M.*t34.*t36.*t40.*u2.*u3+M.*t34.*t36.*t85.*u2.*u3-I2.*I3.*M.*t34.*t36.*u2.*u3.*2.0-I2.*M.*f3.*r.*t8.*t12.*t32-I3.*M.*f3.*r.*t8.*t12.*t34-I2.*M.*f2.*r.*t27.*t34.*t36+I3.*M.*f2.*r.*t27.*t34.*t36+I2.*M.*l.*r.*t8.*t10.*tau1.*2.0+I3.*M.*l.*r.*t8.*t10.*tau1.*2.0-I3.*M.*l.*r.*t8.*t12.*tau2+I2.*M.*t8.*t10.*t27.*t34.*tau3-I3.*M.*t10.*t12.*t34.*t36.*tau2-M.*l.*r.*t27.*t85.*u1.*u2+I2.*t33.*t35.*t36.*t37.*u2.*u3-I3.*t33.*t35.*t36.*t37.*u2.*u3-M.*t34.*t36.*t37.*t40.*u2.*u3+l.*r.*t8.*t10.*t32.*t33.*tau1.*4.0+l.*r.*t8.*t10.*t33.*t34.*tau1.*2.0-l.*r.*t8.*t12.*t32.*t33.*tau2-l.*r.*t8.*t12.*t33.*t34.*tau2+t8.*t10.*t27.*t32.*t33.*t34.*tau3.*3.0-t10.*t12.*t32.*t33.*t34.*t36.*tau2.*3.0+I1.*I2.*M.*l.*r.*t27.*u1.*u2+I2.*I3.*M.*l.*r.*t27.*u1.*u2+I2.*I3.*M.*t34.*t36.*t37.*u2.*u3+I2.*M.*f1.*l.*t8.*t12.*t27.*t34-I3.*M.*f2.*l.*t8.*t10.*t27.*t34.*2.0-I3.*M.*f1.*l.*t8.*t12.*t27.*t34-I2.*M.*f3.*l.*t10.*t12.*t34.*t36.*2.0-I2.*M.*f3.*r.*t8.*t12.*t34.*t36+I3.*M.*f3.*r.*t8.*t12.*t34.*t36+I2.*M.*f2.*r.*t27.*t34.*t36.*t37-I3.*M.*f2.*r.*t27.*t34.*t36.*t37+I2.*g.*r.*t8.*t12.*t27.*t32.*t33-I3.*g.*r.*t8.*t12.*t27.*t32.*t33+I1.*l.*r.*t27.*t32.*t33.*u1.*u2-I2.*l.*r.*t27.*t32.*t33.*u1.*u2+I1.*l.*r.*t27.*t33.*t34.*u1.*u2+I3.*l.*r.*t27.*t32.*t33.*u1.*u2-I2.*l.*r.*t27.*t33.*t34.*u1.*u2+I3.*l.*r.*t27.*t33.*t34.*u1.*u2-M.*l.*r.*t8.*t10.*t40.*u2.*u3.*2.0-M.*l.*r.*t8.*t12.*t40.*u1.*u3+M.*l.*r.*t8.*t10.*t85.*u2.*u3.*2.0+I2.*t8.*t12.*t27.*t32.*t33.*t34.*t41+I2.*t8.*t12.*t27.*t32.*t33.*t34.*t42-I3.*t8.*t12.*t27.*t32.*t33.*t34.*t41-I3.*t8.*t12.*t27.*t32.*t33.*t34.*t42+I1.*t8.*t10.*t27.*t33.*t35.*u1.*u2-I2.*t8.*t10.*t27.*t33.*t35.*u1.*u2+I1.*t10.*t12.*t33.*t35.*t36.*u1.*u3-I3.*t10.*t12.*t33.*t35.*t36.*u1.*u3+I2.*t32.*t33.*t34.*t36.*t37.*u2.*u3.*5.0-I3.*t32.*t33.*t34.*t36.*t37.*u2.*u3.*5.0-M.*t10.*t12.*t34.*t36.*t40.*u1.*u3-M.*t8.*t10.*t27.*t34.*t85.*u1.*u2+l.*r.*t27.*t33.*t34.*t36.*t37.*tau3.*2.0+I2.*M.*f1.*r.*t10.*t12.*t27.*t34.*t36-I3.*M.*f1.*r.*t10.*t12.*t27.*t34.*t36+I2.*g.*l.*t10.*t12.*t27.*t33.*t34.*t36-I3.*g.*l.*t10.*t12.*t27.*t33.*t34.*t36+I1.*l.*r.*t8.*t12.*t32.*t33.*u1.*u3+I2.*l.*r.*t8.*t10.*t32.*t33.*u2.*u3.*4.0+I2.*l.*r.*t8.*t12.*t32.*t33.*u1.*u3-I3.*l.*r.*t8.*t10.*t32.*t33.*u2.*u3.*4.0+I1.*l.*r.*t8.*t12.*t33.*t34.*u1.*u3+I2.*l.*r.*t8.*t10.*t33.*t34.*u2.*u3.*2.0-I3.*l.*r.*t8.*t12.*t32.*t33.*u1.*u3-I3.*l.*r.*t8.*t10.*t33.*t34.*u2.*u3.*2.0+I2.*l.*r.*t27.*t33.*t34.*t36.*u1.*u2-I3.*l.*r.*t27.*t33.*t34.*t36.*u1.*u2+I1.*t8.*t10.*t27.*t32.*t33.*t34.*u1.*u2.*3.0-I2.*t8.*t10.*t27.*t32.*t33.*t34.*u1.*u2.*3.0+I3.*t8.*t10.*t27.*t32.*t33.*t34.*u1.*u2.*2.0+I1.*t10.*t12.*t32.*t33.*t34.*t36.*u1.*u3.*3.0+I2.*t10.*t12.*t32.*t33.*t34.*t36.*u1.*u3.*2.0-I3.*t10.*t12.*t32.*t33.*t34.*t36.*u1.*u3.*3.0+l.*r.*t8.*t10.*t33.*t34.*t36.*t37.*tau1.*2.0-l.*r.*t8.*t12.*t33.*t34.*t36.*t37.*tau2.*2.0+I1.*I3.*M.*l.*r.*t8.*t12.*u1.*u3+I2.*I3.*M.*l.*r.*t8.*t12.*u1.*u3+I1.*I2.*M.*t8.*t10.*t27.*t34.*u1.*u2+I1.*I3.*M.*t10.*t12.*t34.*t36.*u1.*u3+I2.*l.*r.*t10.*t12.*t27.*t33.*t34.*t36.*t41+I2.*l.*r.*t10.*t12.*t27.*t33.*t34.*t36.*t42-I3.*l.*r.*t10.*t12.*t27.*t33.*t34.*t36.*t41-I3.*l.*r.*t10.*t12.*t27.*t33.*t34.*t36.*t42+I2.*l.*r.*t8.*t12.*t33.*t34.*t36.*u1.*u3-I3.*l.*r.*t8.*t12.*t33.*t34.*t36.*u1.*u3+I1.*l.*r.*t27.*t33.*t34.*t36.*t37.*u1.*u2.*2.0-I2.*l.*r.*t27.*t33.*t34.*t36.*t37.*u1.*u2.*3.0+I3.*l.*r.*t27.*t33.*t34.*t36.*t37.*u1.*u2+I1.*l.*r.*t8.*t12.*t33.*t34.*t36.*t37.*u1.*u3.*2.0+I2.*l.*r.*t8.*t10.*t33.*t34.*t36.*t37.*u2.*u3.*2.0-I3.*l.*r.*t8.*t10.*t33.*t34.*t36.*t37.*u2.*u3.*2.0-I3.*l.*r.*t8.*t12.*t33.*t34.*t36.*t37.*u1.*u3.*2.0)-t100.*(I2.*I3.*tau1p1+t33.*t39.*tau1p1-I2.*t40.*u2p1.*u3p1+I3.*t85.*u2p1.*u3p1+t32.*t33.*t34.*tau1p1+I2.*M.*t32.*tau1p1+I3.*M.*t32.*tau1p1+I3.*M.*t34.*tau1p1-I2.*I3.*f2p1.*r.*t5+I2.*M.*t34.*t50.*tau1p1-I3.*M.*t34.*t50.*tau1p1+I2.*t33.*t39.*u2p1.*u3p1-I3.*t33.*t39.*u2p1.*u3p1-M.*t32.*t40.*u2p1.*u3p1-M.*t34.*t40.*u2p1.*u3p1+M.*t32.*t85.*u2p1.*u3p1+t33.*t35.*t50.*t53.*tau1p1+l.*r.*t5.*t32.*t33.*tau3p1+l.*r.*t5.*t33.*t34.*tau3p1+t3.*t5.*t7.*t33.*t35.*tau3p1-t3.*t6.*t33.*t35.*t50.*tau2p1+t32.*t33.*t34.*t50.*t53.*tau1p1.*5.0+I2.*I3.*M.*t34.*u2p1.*u3p1-I2.*I3.*f3p1.*r.*t6.*t7-I3.*M.*f2p1.*r.*t5.*t32-I3.*M.*f2p1.*r.*t5.*t34+I2.*M.*l.*r.*t5.*tau3p1+I3.*M.*t34.*t50.*t53.*tau1p1+I2.*t32.*t33.*t34.*u2p1.*u3p1-I3.*t32.*t33.*t34.*u2p1.*u3p1+M.*t34.*t40.*t50.*u2p1.*u3p1+M.*t34.*t50.*t85.*u2p1.*u3p1-I2.*I3.*M.*t34.*t50.*u2p1.*u3p1.*2.0-I2.*M.*f3p1.*r.*t6.*t7.*t32-I3.*M.*f3p1.*r.*t6.*t7.*t34-I2.*M.*f2p1.*r.*t5.*t34.*t50+I3.*M.*f2p1.*r.*t5.*t34.*t50+I2.*M.*l.*r.*t3.*t7.*tau1p1.*2.0+I3.*M.*l.*r.*t3.*t7.*tau1p1.*2.0-I3.*M.*l.*r.*t6.*t7.*tau2p1+I2.*M.*t3.*t5.*t7.*t34.*tau3p1-I3.*M.*t3.*t6.*t34.*t50.*tau2p1-M.*l.*r.*t5.*t85.*u1p1.*u2p1+I2.*t33.*t35.*t50.*t53.*u2p1.*u3p1-I3.*t33.*t35.*t50.*t53.*u2p1.*u3p1-M.*t34.*t40.*t50.*t53.*u2p1.*u3p1+l.*r.*t3.*t7.*t32.*t33.*tau1p1.*4.0+l.*r.*t3.*t7.*t33.*t34.*tau1p1.*2.0-l.*r.*t6.*t7.*t32.*t33.*tau2p1-l.*r.*t6.*t7.*t33.*t34.*tau2p1+t3.*t5.*t7.*t32.*t33.*t34.*tau3p1.*3.0-t3.*t6.*t32.*t33.*t34.*t50.*tau2p1.*3.0+I1.*I2.*M.*l.*r.*t5.*u1p1.*u2p1+I2.*I3.*M.*l.*r.*t5.*u1p1.*u2p1+I2.*I3.*M.*t34.*t50.*t53.*u2p1.*u3p1+I2.*M.*f1p1.*l.*t5.*t6.*t7.*t34-I3.*M.*f1p1.*l.*t5.*t6.*t7.*t34-I3.*M.*f2p1.*l.*t3.*t5.*t7.*t34.*2.0-I2.*M.*f3p1.*l.*t3.*t6.*t34.*t50.*2.0-I2.*M.*f3p1.*r.*t6.*t7.*t34.*t50+I3.*M.*f3p1.*r.*t6.*t7.*t34.*t50+I2.*M.*f2p1.*r.*t5.*t34.*t50.*t53-I3.*M.*f2p1.*r.*t5.*t34.*t50.*t53+I2.*g.*r.*t5.*t6.*t7.*t32.*t33-I3.*g.*r.*t5.*t6.*t7.*t32.*t33+I1.*l.*r.*t5.*t32.*t33.*u1p1.*u2p1-I2.*l.*r.*t5.*t32.*t33.*u1p1.*u2p1+I1.*l.*r.*t5.*t33.*t34.*u1p1.*u2p1+I3.*l.*r.*t5.*t32.*t33.*u1p1.*u2p1-I2.*l.*r.*t5.*t33.*t34.*u1p1.*u2p1+I3.*l.*r.*t5.*t33.*t34.*u1p1.*u2p1-M.*l.*r.*t6.*t7.*t40.*u1p1.*u3p1-M.*l.*r.*t3.*t7.*t40.*u2p1.*u3p1.*2.0+M.*l.*r.*t3.*t7.*t85.*u2p1.*u3p1.*2.0+I2.*t5.*t6.*t7.*t32.*t33.*t34.*t54+I2.*t5.*t6.*t7.*t32.*t33.*t34.*t55-I3.*t5.*t6.*t7.*t32.*t33.*t34.*t54-I3.*t5.*t6.*t7.*t32.*t33.*t34.*t55+I1.*t3.*t5.*t7.*t33.*t35.*u1p1.*u2p1-I2.*t3.*t5.*t7.*t33.*t35.*u1p1.*u2p1+I1.*t3.*t6.*t33.*t35.*t50.*u1p1.*u3p1-I3.*t3.*t6.*t33.*t35.*t50.*u1p1.*u3p1+I2.*t32.*t33.*t34.*t50.*t53.*u2p1.*u3p1.*5.0-I3.*t32.*t33.*t34.*t50.*t53.*u2p1.*u3p1.*5.0-M.*t3.*t5.*t7.*t34.*t85.*u1p1.*u2p1-M.*t3.*t6.*t34.*t40.*t50.*u1p1.*u3p1+l.*r.*t5.*t33.*t34.*t50.*t53.*tau3p1.*2.0+I2.*M.*f1p1.*r.*t3.*t5.*t6.*t34.*t50-I3.*M.*f1p1.*r.*t3.*t5.*t6.*t34.*t50+I2.*g.*l.*t3.*t5.*t6.*t33.*t34.*t50-I3.*g.*l.*t3.*t5.*t6.*t33.*t34.*t50+I1.*l.*r.*t6.*t7.*t32.*t33.*u1p1.*u3p1+I2.*l.*r.*t6.*t7.*t32.*t33.*u1p1.*u3p1+I1.*l.*r.*t6.*t7.*t33.*t34.*u1p1.*u3p1-I3.*l.*r.*t6.*t7.*t32.*t33.*u1p1.*u3p1+I2.*l.*r.*t3.*t7.*t32.*t33.*u2p1.*u3p1.*4.0-I3.*l.*r.*t3.*t7.*t32.*t33.*u2p1.*u3p1.*4.0+I2.*l.*r.*t3.*t7.*t33.*t34.*u2p1.*u3p1.*2.0-I3.*l.*r.*t3.*t7.*t33.*t34.*u2p1.*u3p1.*2.0+I2.*l.*r.*t5.*t33.*t34.*t50.*u1p1.*u2p1-I3.*l.*r.*t5.*t33.*t34.*t50.*u1p1.*u2p1+I1.*t3.*t5.*t7.*t32.*t33.*t34.*u1p1.*u2p1.*3.0-I2.*t3.*t5.*t7.*t32.*t33.*t34.*u1p1.*u2p1.*3.0+I3.*t3.*t5.*t7.*t32.*t33.*t34.*u1p1.*u2p1.*2.0+I1.*t3.*t6.*t32.*t33.*t34.*t50.*u1p1.*u3p1.*3.0+I2.*t3.*t6.*t32.*t33.*t34.*t50.*u1p1.*u3p1.*2.0-I3.*t3.*t6.*t32.*t33.*t34.*t50.*u1p1.*u3p1.*3.0+l.*r.*t3.*t7.*t33.*t34.*t50.*t53.*tau1p1.*2.0-l.*r.*t6.*t7.*t33.*t34.*t50.*t53.*tau2p1.*2.0+I1.*I3.*M.*l.*r.*t6.*t7.*u1p1.*u3p1+I2.*I3.*M.*l.*r.*t6.*t7.*u1p1.*u3p1+I1.*I2.*M.*t3.*t5.*t7.*t34.*u1p1.*u2p1+I1.*I3.*M.*t3.*t6.*t34.*t50.*u1p1.*u3p1+I2.*l.*r.*t3.*t5.*t6.*t33.*t34.*t50.*t54+I2.*l.*r.*t3.*t5.*t6.*t33.*t34.*t50.*t55-I3.*l.*r.*t3.*t5.*t6.*t33.*t34.*t50.*t54-I3.*l.*r.*t3.*t5.*t6.*t33.*t34.*t50.*t55+I2.*l.*r.*t6.*t7.*t33.*t34.*t50.*u1p1.*u3p1-I3.*l.*r.*t6.*t7.*t33.*t34.*t50.*u1p1.*u3p1+I1.*l.*r.*t5.*t33.*t34.*t50.*t53.*u1p1.*u2p1.*2.0-I2.*l.*r.*t5.*t33.*t34.*t50.*t53.*u1p1.*u2p1.*3.0+I3.*l.*r.*t5.*t33.*t34.*t50.*t53.*u1p1.*u2p1+I1.*l.*r.*t6.*t7.*t33.*t34.*t50.*t53.*u1p1.*u3p1.*2.0-I3.*l.*r.*t6.*t7.*t33.*t34.*t50.*t53.*u1p1.*u3p1.*2.0+I2.*l.*r.*t3.*t7.*t33.*t34.*t50.*t53.*u2p1.*u3p1.*2.0-I3.*l.*r.*t3.*t7.*t33.*t34.*t50.*t53.*u2p1.*u3p1.*2.0)).*(1.0./8.0))-t60.*cos(t68).*(u3.*(1.0./2.0)+u3p1.*(1.0./2.0)+dt.*(t84.*(I1.*I2.*tau3+t33.*t35.*tau3+I2.*t38.*u1.*u2-I1.*t85.*u1.*u2+t32.*t33.*t34.*tau3-t33.*t35.*t36.*tau3+I1.*I2.*f2.*l+I1.*M.*t32.*tau3+I1.*M.*t34.*tau3+I2.*M.*t34.*tau3+I1.*M.*f2.*l.*t32+I1.*M.*f2.*l.*t34-I1.*M.*t34.*t36.*tau3+I1.*t33.*t35.*u1.*u2-I2.*t33.*t35.*u1.*u2-I1.*t33.*t39.*u1.*u2+M.*t32.*t38.*u1.*u2+M.*t34.*t38.*u1.*u2-M.*t34.*t85.*u1.*u2-t32.*t33.*t34.*t36.*tau3+l.*r.*t27.*t32.*t33.*tau1+l.*r.*t27.*t33.*t34.*tau1+t8.*t10.*t27.*t33.*t35.*tau1-t8.*t12.*t27.*t33.*t35.*tau2-I1.*I2.*M.*t32.*u1.*u2.*2.0-I1.*M.*f2.*l.*t34.*t36+I2.*M.*f2.*l.*t34.*t36+I1.*I2.*f2.*r.*t8.*t10+I1.*I2.*f1.*r.*t8.*t12+I2.*M.*l.*r.*t27.*tau1+I1.*M.*t34.*t36.*t37.*tau3-I2.*M.*t34.*t36.*t37.*tau3-I2.*t32.*t33.*t34.*u1.*u2-I1.*t33.*t35.*t36.*u1.*u2+I2.*t33.*t35.*t36.*u1.*u2-M.*t34.*t36.*t38.*u1.*u2+I1.*I2.*M.*t34.*t36.*u1.*u2+I1.*M.*f2.*l.*t34.*t36.*t37.*3.0-I2.*M.*f2.*l.*t34.*t36.*t37+I1.*M.*f2.*r.*t8.*t10.*t32.*3.0+I1.*M.*f1.*r.*t8.*t12.*t32+I1.*M.*f2.*r.*t8.*t10.*t34+I2.*M.*f1.*r.*t8.*t12.*t34+I1.*M.*l.*r.*t8.*t10.*tau3.*2.0+I2.*M.*t8.*t10.*t27.*t34.*tau1-I1.*M.*t8.*t12.*t27.*t34.*tau2+I1.*g.*l.*t8.*t12.*t32.*t33+I2.*g.*l.*t8.*t12.*t33.*t34+M.*l.*r.*t27.*t85.*u2.*u3+M.*t34.*t36.*t37.*t38.*u1.*u2+M.*t34.*t36.*t37.*t85.*u1.*u2+l.*r.*t8.*t10.*t33.*t34.*tau3.*2.0+t8.*t10.*t27.*t32.*t33.*t34.*tau1.*3.0-t8.*t12.*t27.*t32.*t33.*t34.*tau2+I1.*I2.*M.*g.*l.*t8.*t12+I1.*I2.*M.*l.*r.*t8.*t12.*t41+I1.*I2.*M.*l.*r.*t8.*t12.*t42-I2.*I3.*M.*l.*r.*t27.*u2.*u3-I1.*I2.*M.*t34.*t36.*t37.*u1.*u2.*2.0+I1.*M.*f3.*l.*t8.*t12.*t27.*t34-I2.*M.*f3.*l.*t8.*t12.*t27.*t34+I1.*M.*f1.*l.*t10.*t12.*t34.*t36.*2.0-I1.*M.*f2.*r.*t8.*t10.*t34.*t36+I2.*M.*f2.*r.*t8.*t10.*t34.*t36+I1.*g.*r.*t10.*t12.*t32.*t33.*t36.*2.0+I1.*l.*r.*t8.*t12.*t32.*t33.*t41+I1.*l.*r.*t8.*t12.*t32.*t33.*t42+I2.*l.*r.*t8.*t12.*t33.*t34.*t41+I2.*l.*r.*t8.*t12.*t33.*t34.*t42+I2.*l.*r.*t27.*t32.*t33.*u2.*u3-I3.*l.*r.*t27.*t32.*t33.*u2.*u3+I2.*l.*r.*t27.*t33.*t34.*u2.*u3-I3.*l.*r.*t27.*t33.*t34.*u2.*u3+M.*l.*r.*t8.*t10.*t38.*u1.*u2.*2.0+I1.*t10.*t12.*t32.*t33.*t34.*t36.*t41.*2.0+I1.*t10.*t12.*t32.*t33.*t34.*t36.*t42.*2.0+I1.*t8.*t12.*t27.*t33.*t35.*u1.*u3+I2.*t8.*t10.*t27.*t33.*t35.*u2.*u3-I3.*t8.*t10.*t27.*t33.*t35.*u2.*u3-I3.*t8.*t12.*t27.*t33.*t35.*u1.*u3-I1.*t32.*t33.*t34.*t36.*t37.*u1.*u2.*3.0+I2.*t32.*t33.*t34.*t36.*t37.*u1.*u2+M.*t8.*t12.*t27.*t34.*t38.*u1.*u3+M.*t8.*t10.*t27.*t34.*t85.*u2.*u3-l.*r.*t8.*t10.*t33.*t34.*t36.*tau3.*2.0+l.*r.*t27.*t33.*t34.*t36.*t37.*tau1.*2.0+I1.*M.*f3.*r.*t10.*t12.*t27.*t34.*t36-I2.*M.*f3.*r.*t10.*t12.*t27.*t34.*t36+I1.*M.*f2.*r.*t8.*t10.*t34.*t36.*t37+I1.*M.*f1.*r.*t8.*t12.*t34.*t36.*t37-I2.*M.*f2.*r.*t8.*t10.*t34.*t36.*t37-I2.*M.*f1.*r.*t8.*t12.*t34.*t36.*t37+I1.*g.*l.*t8.*t12.*t33.*t34.*t36.*t37-I2.*g.*l.*t8.*t12.*t33.*t34.*t36.*t37-I1.*l.*r.*t8.*t10.*t32.*t33.*u1.*u2.*3.0+I1.*l.*r.*t8.*t10.*t33.*t34.*u1.*u2-I2.*l.*r.*t8.*t10.*t33.*t34.*u1.*u2.*2.0+I2.*t8.*t10.*t27.*t32.*t33.*t34.*u2.*u3.*3.0+I2.*t8.*t12.*t27.*t32.*t33.*t34.*u1.*u3-I3.*t8.*t10.*t27.*t32.*t33.*t34.*u2.*u3.*3.0-I3.*t8.*t12.*t27.*t32.*t33.*t34.*u1.*u3-l.*r.*t10.*t12.*t27.*t33.*t34.*t36.*tau2.*2.0-I1.*I2.*M.*l.*r.*t8.*t10.*u1.*u2.*3.0-I1.*I3.*M.*t8.*t12.*t27.*t34.*u1.*u3-I2.*I3.*M.*t8.*t10.*t27.*t34.*u2.*u3+I1.*l.*r.*t8.*t12.*t33.*t34.*t36.*t37.*t41+I1.*l.*r.*t8.*t12.*t33.*t34.*t36.*t37.*t42-I2.*l.*r.*t8.*t12.*t33.*t34.*t36.*t37.*t41-I2.*l.*r.*t8.*t12.*t33.*t34.*t36.*t37.*t42-I1.*l.*r.*t8.*t10.*t33.*t34.*t36.*u1.*u2+I2.*l.*r.*t8.*t10.*t33.*t34.*t36.*u1.*u2+I2.*l.*r.*t27.*t33.*t34.*t36.*t37.*u2.*u3.*2.0-I3.*l.*r.*t27.*t33.*t34.*t36.*t37.*u2.*u3.*2.0+I1.*l.*r.*t10.*t12.*t27.*t33.*t34.*t36.*u1.*u3+I2.*l.*r.*t10.*t12.*t27.*t33.*t34.*t36.*u1.*u3-I3.*l.*r.*t10.*t12.*t27.*t33.*t34.*t36.*u1.*u3.*2.0-I1.*l.*r.*t8.*t10.*t33.*t34.*t36.*t37.*u1.*u2+I2.*l.*r.*t8.*t10.*t33.*t34.*t36.*t37.*u1.*u2)-t100.*(I1.*I2.*tau3p1+t33.*t35.*tau3p1+I2.*t38.*u1p1.*u2p1-I1.*t85.*u1p1.*u2p1+t32.*t33.*t34.*tau3p1-t33.*t35.*t50.*tau3p1+I1.*I2.*f2p1.*l+I1.*M.*t32.*tau3p1+I1.*M.*t34.*tau3p1+I2.*M.*t34.*tau3p1+I1.*M.*f2p1.*l.*t32+I1.*M.*f2p1.*l.*t34-I1.*M.*t34.*t50.*tau3p1+I1.*t33.*t35.*u1p1.*u2p1-I2.*t33.*t35.*u1p1.*u2p1-I1.*t33.*t39.*u1p1.*u2p1+M.*t32.*t38.*u1p1.*u2p1+M.*t34.*t38.*u1p1.*u2p1-M.*t34.*t85.*u1p1.*u2p1-t32.*t33.*t34.*t50.*tau3p1+l.*r.*t5.*t32.*t33.*tau1p1+l.*r.*t5.*t33.*t34.*tau1p1+t3.*t5.*t7.*t33.*t35.*tau1p1-t5.*t6.*t7.*t33.*t35.*tau2p1-I1.*I2.*M.*t32.*u1p1.*u2p1.*2.0-I1.*M.*f2p1.*l.*t34.*t50+I2.*M.*f2p1.*l.*t34.*t50+I1.*I2.*f1p1.*r.*t6.*t7+I1.*I2.*f2p1.*r.*t3.*t7+I2.*M.*l.*r.*t5.*tau1p1+I1.*M.*t34.*t50.*t53.*tau3p1-I2.*M.*t34.*t50.*t53.*tau3p1-I2.*t32.*t33.*t34.*u1p1.*u2p1-I1.*t33.*t35.*t50.*u1p1.*u2p1+I2.*t33.*t35.*t50.*u1p1.*u2p1-M.*t34.*t38.*t50.*u1p1.*u2p1+I1.*I2.*M.*t34.*t50.*u1p1.*u2p1+I1.*M.*f2p1.*l.*t34.*t50.*t53.*3.0-I2.*M.*f2p1.*l.*t34.*t50.*t53+I1.*M.*f1p1.*r.*t6.*t7.*t32+I2.*M.*f1p1.*r.*t6.*t7.*t34+I1.*M.*f2p1.*r.*t3.*t7.*t32.*3.0+I1.*M.*f2p1.*r.*t3.*t7.*t34+I1.*M.*l.*r.*t3.*t7.*tau3p1.*2.0+I2.*M.*t3.*t5.*t7.*t34.*tau1p1-I1.*M.*t5.*t6.*t7.*t34.*tau2p1+I1.*g.*l.*t6.*t7.*t32.*t33+I2.*g.*l.*t6.*t7.*t33.*t34+M.*l.*r.*t5.*t85.*u2p1.*u3p1+M.*t34.*t38.*t50.*t53.*u1p1.*u2p1+M.*t34.*t50.*t53.*t85.*u1p1.*u2p1+l.*r.*t3.*t7.*t33.*t34.*tau3p1.*2.0+t3.*t5.*t7.*t32.*t33.*t34.*tau1p1.*3.0-t5.*t6.*t7.*t32.*t33.*t34.*tau2p1+I1.*I2.*M.*g.*l.*t6.*t7+I1.*I2.*M.*l.*r.*t6.*t7.*t54+I1.*I2.*M.*l.*r.*t6.*t7.*t55-I2.*I3.*M.*l.*r.*t5.*u2p1.*u3p1-I1.*I2.*M.*t34.*t50.*t53.*u1p1.*u2p1.*2.0+I1.*M.*f3p1.*l.*t5.*t6.*t7.*t34-I2.*M.*f3p1.*l.*t5.*t6.*t7.*t34+I1.*M.*f1p1.*l.*t3.*t6.*t34.*t50.*2.0-I1.*M.*f2p1.*r.*t3.*t7.*t34.*t50+I2.*M.*f2p1.*r.*t3.*t7.*t34.*t50+I1.*g.*r.*t3.*t6.*t32.*t33.*t50.*2.0+I1.*l.*r.*t6.*t7.*t32.*t33.*t54+I1.*l.*r.*t6.*t7.*t32.*t33.*t55+I2.*l.*r.*t6.*t7.*t33.*t34.*t54+I2.*l.*r.*t6.*t7.*t33.*t34.*t55+I2.*l.*r.*t5.*t32.*t33.*u2p1.*u3p1-I3.*l.*r.*t5.*t32.*t33.*u2p1.*u3p1+I2.*l.*r.*t5.*t33.*t34.*u2p1.*u3p1-I3.*l.*r.*t5.*t33.*t34.*u2p1.*u3p1+M.*l.*r.*t3.*t7.*t38.*u1p1.*u2p1.*2.0+I1.*t3.*t6.*t32.*t33.*t34.*t50.*t54.*2.0+I1.*t3.*t6.*t32.*t33.*t34.*t50.*t55.*2.0+I1.*t5.*t6.*t7.*t33.*t35.*u1p1.*u3p1-I3.*t5.*t6.*t7.*t33.*t35.*u1p1.*u3p1+I2.*t3.*t5.*t7.*t33.*t35.*u2p1.*u3p1-I3.*t3.*t5.*t7.*t33.*t35.*u2p1.*u3p1-I1.*t32.*t33.*t34.*t50.*t53.*u1p1.*u2p1.*3.0+I2.*t32.*t33.*t34.*t50.*t53.*u1p1.*u2p1+M.*t5.*t6.*t7.*t34.*t38.*u1p1.*u3p1+M.*t3.*t5.*t7.*t34.*t85.*u2p1.*u3p1-l.*r.*t3.*t7.*t33.*t34.*t50.*tau3p1.*2.0+l.*r.*t5.*t33.*t34.*t50.*t53.*tau1p1.*2.0+I1.*M.*f3p1.*r.*t3.*t5.*t6.*t34.*t50-I2.*M.*f3p1.*r.*t3.*t5.*t6.*t34.*t50+I1.*M.*f1p1.*r.*t6.*t7.*t34.*t50.*t53-I2.*M.*f1p1.*r.*t6.*t7.*t34.*t50.*t53+I1.*M.*f2p1.*r.*t3.*t7.*t34.*t50.*t53-I2.*M.*f2p1.*r.*t3.*t7.*t34.*t50.*t53+I1.*g.*l.*t6.*t7.*t33.*t34.*t50.*t53-I2.*g.*l.*t6.*t7.*t33.*t34.*t50.*t53-I1.*l.*r.*t3.*t7.*t32.*t33.*u1p1.*u2p1.*3.0+I1.*l.*r.*t3.*t7.*t33.*t34.*u1p1.*u2p1-I2.*l.*r.*t3.*t7.*t33.*t34.*u1p1.*u2p1.*2.0+I2.*t5.*t6.*t7.*t32.*t33.*t34.*u1p1.*u3p1-I3.*t5.*t6.*t7.*t32.*t33.*t34.*u1p1.*u3p1+I2.*t3.*t5.*t7.*t32.*t33.*t34.*u2p1.*u3p1.*3.0-I3.*t3.*t5.*t7.*t32.*t33.*t34.*u2p1.*u3p1.*3.0-l.*r.*t3.*t5.*t6.*t33.*t34.*t50.*tau2p1.*2.0-I1.*I2.*M.*l.*r.*t3.*t7.*u1p1.*u2p1.*3.0-I1.*I3.*M.*t5.*t6.*t7.*t34.*u1p1.*u3p1-I2.*I3.*M.*t3.*t5.*t7.*t34.*u2p1.*u3p1+I1.*l.*r.*t6.*t7.*t33.*t34.*t50.*t53.*t54+I1.*l.*r.*t6.*t7.*t33.*t34.*t50.*t53.*t55-I2.*l.*r.*t6.*t7.*t33.*t34.*t50.*t53.*t54-I2.*l.*r.*t6.*t7.*t33.*t34.*t50.*t53.*t55-I1.*l.*r.*t3.*t7.*t33.*t34.*t50.*u1p1.*u2p1+I2.*l.*r.*t3.*t7.*t33.*t34.*t50.*u1p1.*u2p1+I2.*l.*r.*t5.*t33.*t34.*t50.*t53.*u2p1.*u3p1.*2.0-I3.*l.*r.*t5.*t33.*t34.*t50.*t53.*u2p1.*u3p1.*2.0+I1.*l.*r.*t3.*t5.*t6.*t33.*t34.*t50.*u1p1.*u3p1+I2.*l.*r.*t3.*t5.*t6.*t33.*t34.*t50.*u1p1.*u3p1-I3.*l.*r.*t3.*t5.*t6.*t33.*t34.*t50.*u1p1.*u3p1.*2.0-I1.*l.*r.*t3.*t7.*t33.*t34.*t50.*t53.*u1p1.*u2p1+I2.*l.*r.*t3.*t7.*t33.*t34.*t50.*t53.*u1p1.*u2p1)).*(1.0./8.0)).*(1.0./2.0))-r.*(u1p1.*(t2.*t6+t3.*t4.*t5)+u2p1.*(t2.*t3-t4.*t5.*t6)-t4.*t7.*u3p1).*(1.0./4.0);
