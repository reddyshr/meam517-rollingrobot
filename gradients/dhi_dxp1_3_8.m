function out1 = dhi_dxp1_3_8(q1,q2,q3,q4,q5,u1,u2,u3,tau1,tau2,tau3,f1,f2,f3,q1p1,q2p1,q3p1,q4p1,q5p1,u1p1,u2p1,u3p1,tau1p1,tau2p1,tau3p1,f1p1,f2p1,f3p1,dt,M,l,r,g,I1,I2,I3)
%DHI_DXP1_3_8
%    OUT1 = DHI_DXP1_3_8(Q1,Q2,Q3,Q4,Q5,U1,U2,U3,TAU1,TAU2,TAU3,F1,F2,F3,Q1P1,Q2P1,Q3P1,Q4P1,Q5P1,U1P1,U2P1,U3P1,TAU1P1,TAU2P1,TAU3P1,F1P1,F2P1,F3P1,DT,M,L,R,G,I1,I2,I3)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    23-Nov-2018 20:52:51

t2 = r.^2;
t3 = sin(q2p1);
t4 = M.^2;
t5 = l.^2;
t6 = I2.^2;
t7 = cos(q2p1);
t8 = cos(q3p1);
t9 = t2.^2;
t10 = sin(q3p1);
t11 = t7.^2;
t12 = t8.^2;
t13 = q2.*(1.0./2.0);
t14 = q2p1.*(1.0./2.0);
t15 = cos(q3);
t16 = t15.*u2;
t17 = sin(q3);
t18 = t17.*u1;
t19 = t16+t18-t10.*u1p1-t8.*u2p1;
t20 = dt.*t19.*(1.0./8.0);
t21 = t13+t14+t20;
t22 = t5.^2;
t23 = I1.*t4.*t22;
t24 = I3.*t4.*t9;
t25 = I1.*I2.*I3;
t26 = I1.*I2.*M.*t5;
t27 = I1.*I3.*M.*t5;
t28 = I1.*I3.*M.*t2;
t29 = I2.*I3.*M.*t2;
t30 = cos(q2);
t31 = t30.^2;
t32 = I1.*t2.*t4.*t5;
t33 = I3.*t2.*t4.*t5;
t34 = t15.^2;
t35 = I3.^2;
t36 = sin(q2);
t37 = u2.^2;
t38 = u3.^2;
t39 = I2.*t4.*t9.*t11;
t40 = I2.*t2.*t4.*t5.*t11;
t41 = I1.*I2.*M.*t2.*t11;
t42 = I1.*t4.*t9.*t11.*t12;
t43 = I1.*l.*r.*t2.*t4.*t7.*t8.*2.0;
t44 = I1.*l.*r.*t4.*t5.*t7.*t8.*4.0;
t45 = I3.*l.*r.*t2.*t4.*t7.*t8.*2.0;
t46 = I1.*t2.*t4.*t5.*t11.*t12.*5.0;
t47 = I1.*I3.*M.*t2.*t11.*t12;
t48 = I2.*l.*r.*t2.*t4.*t7.*t8.*t11.*2.0;
t49 = I1.*l.*r.*t2.*t4.*t7.*t8.*t11.*t12.*2.0;
t50 = I1.*I2.*M.*l.*r.*t7.*t8.*2.0;
t51 = I1.*I3.*M.*l.*r.*t7.*t8.*2.0;
t84 = I3.*t4.*t9.*t11;
t85 = I3.*t2.*t4.*t5.*t11;
t86 = I1.*I3.*M.*t2.*t11;
t87 = I2.*t4.*t9.*t11.*t12;
t88 = I2.*t2.*t4.*t5.*t11.*t12;
t89 = I2.*I3.*M.*t2.*t11.*t12;
t90 = I3.*l.*r.*t2.*t4.*t7.*t8.*t11.*2.0;
t91 = I2.*l.*r.*t2.*t4.*t7.*t8.*t11.*t12.*2.0;
t52 = t23+t24+t25+t26+t27+t28+t29+t32+t33+t39+t40+t41+t42+t43+t44+t45+t46+t47+t48+t49+t50+t51-t84-t85-t86-t87-t88-t89-t90-t91;
t53 = 1.0./t52;
t54 = u2p1.^2;
t55 = u3p1.^2;
t56 = q3.*(1.0./2.0);
t57 = q3p1.*(1.0./2.0);
t58 = 1.0./t30;
t59 = t15.*u1;
t92 = t17.*u2;
t60 = t59-t92;
t61 = 1.0./t7;
t62 = t8.*u1p1;
t94 = t10.*u2p1;
t63 = t62-t94;
t64 = t3.*t61.*t63;
t93 = t36.*t58.*t60;
t65 = t64-t93+u3-u3p1;
t66 = dt.*t65.*(1.0./8.0);
t67 = t56+t57+t66;
t68 = I2.*t4.*t9.*t31;
t69 = I2.*t2.*t4.*t5.*t31;
t70 = I1.*I2.*M.*t2.*t31;
t71 = I1.*t4.*t9.*t31.*t34;
t72 = I1.*l.*r.*t2.*t4.*t15.*t30.*2.0;
t73 = I1.*l.*r.*t4.*t5.*t15.*t30.*4.0;
t74 = I3.*l.*r.*t2.*t4.*t15.*t30.*2.0;
t75 = I1.*t2.*t4.*t5.*t31.*t34.*5.0;
t76 = I1.*I3.*M.*t2.*t31.*t34;
t77 = I2.*l.*r.*t2.*t4.*t15.*t30.*t31.*2.0;
t78 = I1.*l.*r.*t2.*t4.*t15.*t30.*t31.*t34.*2.0;
t79 = I1.*I2.*M.*l.*r.*t15.*t30.*2.0;
t80 = I1.*I3.*M.*l.*r.*t15.*t30.*2.0;
t81 = t23+t24+t25+t26+t27+t28+t29+t32+t33+t68+t69+t70+t71+t72+t73+t74+t75+t76+t77+t78+t79+t80-I3.*t4.*t9.*t31-I1.*I3.*M.*t2.*t31-I3.*t2.*t4.*t5.*t31-I2.*t4.*t9.*t31.*t34-I2.*I3.*M.*t2.*t31.*t34-I2.*t2.*t4.*t5.*t31.*t34-I3.*l.*r.*t2.*t4.*t15.*t30.*t31.*2.0-I2.*l.*r.*t2.*t4.*t15.*t30.*t31.*t34.*2.0;
t82 = 1.0./t81;
t83 = I1.^2;
t95 = sin(t67);
t96 = cos(t67);
out1 = dt.*t53.*(M.*l.*r.*t3.*t6.*u2p1-I2.*I3.*M.*l.*r.*t3.*u2p1+I2.*l.*r.*t2.*t3.*t4.*u2p1-I3.*l.*r.*t2.*t3.*t4.*u2p1+I2.*l.*r.*t3.*t4.*t5.*u2p1-I3.*l.*r.*t3.*t4.*t5.*u2p1+I1.*t3.*t4.*t7.*t9.*t10.*u1p1-I3.*t3.*t4.*t7.*t9.*t10.*u1p1+I2.*t3.*t4.*t7.*t8.*t9.*u2p1-I3.*t3.*t4.*t7.*t8.*t9.*u2p1+M.*t2.*t3.*t6.*t7.*t8.*u2p1+M.*t2.*t3.*t7.*t10.*t83.*u1p1+I1.*I2.*M.*l.*r.*t7.*t10.*u3p1.*2.0-I1.*I3.*M.*t2.*t3.*t7.*t10.*u1p1-I2.*I3.*M.*t2.*t3.*t7.*t8.*u2p1+I2.*l.*r.*t2.*t4.*t7.*t10.*u3p1.*2.0+I1.*l.*r.*t4.*t5.*t7.*t10.*u3p1.*2.0+I2.*t2.*t3.*t4.*t5.*t7.*t10.*u1p1-I3.*t2.*t3.*t4.*t5.*t7.*t10.*u1p1+I2.*t2.*t3.*t4.*t5.*t7.*t8.*u2p1.*3.0-I3.*t2.*t3.*t4.*t5.*t7.*t8.*u2p1.*3.0+I1.*t2.*t4.*t5.*t8.*t10.*t11.*u3p1.*4.0+I2.*l.*r.*t2.*t3.*t4.*t11.*t12.*u2p1.*2.0-I3.*l.*r.*t2.*t3.*t4.*t11.*t12.*u2p1.*2.0+I1.*l.*r.*t2.*t3.*t4.*t8.*t10.*t11.*u1p1+I2.*l.*r.*t2.*t3.*t4.*t8.*t10.*t11.*u1p1-I3.*l.*r.*t2.*t3.*t4.*t8.*t10.*t11.*u1p1.*2.0+I1.*l.*r.*t2.*t4.*t7.*t10.*t11.*t12.*u3p1.*2.0-I2.*l.*r.*t2.*t4.*t7.*t10.*t11.*t12.*u3p1.*2.0).*(1.0./8.0)+(sin(t21).*(dt.*t95.*(u1.*(1.0./2.0)+u1p1.*(1.0./2.0)+dt.*(t82.*(I2.*I3.*tau1+t4.*t22.*tau1+I3.*t6.*u2.*u3-I2.*t35.*u2.*u3+t2.*t4.*t5.*tau1+I3.*M.*t2.*tau1+I2.*M.*t5.*tau1+I3.*M.*t5.*tau1-I2.*I3.*f2.*r.*t36+I2.*M.*t2.*t31.*tau1-I3.*M.*t2.*t31.*tau1+I2.*t4.*t22.*u2.*u3-I3.*t4.*t22.*u2.*u3+M.*t5.*t6.*u2.*u3-M.*t2.*t35.*u2.*u3-M.*t5.*t35.*u2.*u3+t4.*t9.*t31.*t34.*tau1+l.*r.*t2.*t4.*t36.*tau3+l.*r.*t4.*t5.*t36.*tau3+t2.*t4.*t5.*t31.*t34.*tau1.*5.0-t4.*t9.*t15.*t17.*t31.*tau2+t4.*t9.*t15.*t30.*t36.*tau3+I2.*I3.*M.*t2.*u2.*u3-I2.*I3.*f3.*r.*t17.*t30-I3.*M.*f2.*r.*t2.*t36-I3.*M.*f2.*r.*t5.*t36+I2.*M.*l.*r.*t36.*tau3+I3.*M.*t2.*t31.*t34.*tau1+I2.*t2.*t4.*t5.*u2.*u3-I3.*t2.*t4.*t5.*u2.*u3+M.*t2.*t6.*t31.*u2.*u3+M.*t2.*t31.*t35.*u2.*u3-I2.*I3.*M.*t2.*t31.*u2.*u3.*2.0-I3.*M.*f3.*r.*t2.*t17.*t30-I2.*M.*f3.*r.*t5.*t17.*t30-I2.*M.*f2.*r.*t2.*t31.*t36+I3.*M.*f2.*r.*t2.*t31.*t36+I2.*M.*l.*r.*t15.*t30.*tau1.*2.0+I3.*M.*l.*r.*t15.*t30.*tau1.*2.0-I3.*M.*l.*r.*t17.*t30.*tau2-I3.*M.*t2.*t15.*t17.*t31.*tau2+I2.*M.*t2.*t15.*t30.*t36.*tau3-M.*l.*r.*t6.*t36.*u1.*u2+I2.*t4.*t9.*t31.*t34.*u2.*u3-I3.*t4.*t9.*t31.*t34.*u2.*u3-M.*t2.*t31.*t34.*t35.*u2.*u3+l.*r.*t2.*t4.*t15.*t30.*tau1.*2.0-l.*r.*t2.*t4.*t17.*t30.*tau2+l.*r.*t4.*t5.*t15.*t30.*tau1.*4.0-l.*r.*t4.*t5.*t17.*t30.*tau2-t2.*t4.*t5.*t15.*t17.*t31.*tau2.*3.0+t2.*t4.*t5.*t15.*t30.*t36.*tau3.*3.0+I1.*I2.*M.*l.*r.*t36.*u1.*u2+I2.*I3.*M.*l.*r.*t36.*u1.*u2+I2.*I3.*M.*t2.*t31.*t34.*u2.*u3-I2.*M.*f3.*l.*t2.*t15.*t17.*t31.*2.0+I2.*M.*f1.*l.*t2.*t17.*t30.*t36-I3.*M.*f2.*l.*t2.*t15.*t30.*t36.*2.0-I3.*M.*f1.*l.*t2.*t17.*t30.*t36-I2.*M.*f3.*r.*t2.*t17.*t30.*t31+I3.*M.*f3.*r.*t2.*t17.*t30.*t31+I2.*M.*f2.*r.*t2.*t31.*t34.*t36-I3.*M.*f2.*r.*t2.*t31.*t34.*t36+I2.*g.*r.*t4.*t5.*t17.*t30.*t36-I3.*g.*r.*t4.*t5.*t17.*t30.*t36+I1.*l.*r.*t2.*t4.*t36.*u1.*u2-I2.*l.*r.*t2.*t4.*t36.*u1.*u2+I3.*l.*r.*t2.*t4.*t36.*u1.*u2+I1.*l.*r.*t4.*t5.*t36.*u1.*u2-I2.*l.*r.*t4.*t5.*t36.*u1.*u2+I3.*l.*r.*t4.*t5.*t36.*u1.*u2+M.*l.*r.*t6.*t15.*t30.*u2.*u3.*2.0-M.*l.*r.*t15.*t30.*t35.*u2.*u3.*2.0-M.*l.*r.*t17.*t30.*t35.*u1.*u3+I2.*t2.*t4.*t5.*t17.*t30.*t36.*t37+I2.*t2.*t4.*t5.*t17.*t30.*t36.*t38-I3.*t2.*t4.*t5.*t17.*t30.*t36.*t37-I3.*t2.*t4.*t5.*t17.*t30.*t36.*t38+I1.*t4.*t9.*t15.*t17.*t31.*u1.*u3+I2.*t2.*t4.*t5.*t31.*t34.*u2.*u3.*5.0-I3.*t4.*t9.*t15.*t17.*t31.*u1.*u3-I3.*t2.*t4.*t5.*t31.*t34.*u2.*u3.*5.0+I1.*t4.*t9.*t15.*t30.*t36.*u1.*u2-I2.*t4.*t9.*t15.*t30.*t36.*u1.*u2-M.*t2.*t6.*t15.*t30.*t36.*u1.*u2-M.*t2.*t15.*t17.*t31.*t35.*u1.*u3+l.*r.*t2.*t4.*t31.*t34.*t36.*tau3.*2.0+I2.*M.*f1.*r.*t2.*t15.*t17.*t31.*t36-I3.*M.*f1.*r.*t2.*t15.*t17.*t31.*t36+I2.*g.*l.*t2.*t4.*t15.*t17.*t31.*t36-I3.*g.*l.*t2.*t4.*t15.*t17.*t31.*t36+I1.*l.*r.*t2.*t4.*t17.*t30.*u1.*u3+I2.*l.*r.*t2.*t4.*t15.*t30.*u2.*u3.*2.0-I3.*l.*r.*t2.*t4.*t15.*t30.*u2.*u3.*2.0+I1.*l.*r.*t4.*t5.*t17.*t30.*u1.*u3+I2.*l.*r.*t4.*t5.*t15.*t30.*u2.*u3.*4.0+I2.*l.*r.*t4.*t5.*t17.*t30.*u1.*u3-I3.*l.*r.*t4.*t5.*t15.*t30.*u2.*u3.*4.0-I3.*l.*r.*t4.*t5.*t17.*t30.*u1.*u3+I2.*l.*r.*t2.*t4.*t31.*t36.*u1.*u2-I3.*l.*r.*t2.*t4.*t31.*t36.*u1.*u2+I1.*t2.*t4.*t5.*t15.*t17.*t31.*u1.*u3.*3.0+I2.*t2.*t4.*t5.*t15.*t17.*t31.*u1.*u3.*2.0-I3.*t2.*t4.*t5.*t15.*t17.*t31.*u1.*u3.*3.0+I1.*t2.*t4.*t5.*t15.*t30.*t36.*u1.*u2.*3.0-I2.*t2.*t4.*t5.*t15.*t30.*t36.*u1.*u2.*3.0+I3.*t2.*t4.*t5.*t15.*t30.*t36.*u1.*u2.*2.0+l.*r.*t2.*t4.*t15.*t30.*t31.*t34.*tau1.*2.0-l.*r.*t2.*t4.*t17.*t30.*t31.*t34.*tau2.*2.0+I1.*I3.*M.*l.*r.*t17.*t30.*u1.*u3+I2.*I3.*M.*l.*r.*t17.*t30.*u1.*u3+I1.*I3.*M.*t2.*t15.*t17.*t31.*u1.*u3+I1.*I2.*M.*t2.*t15.*t30.*t36.*u1.*u2+I2.*l.*r.*t2.*t4.*t15.*t17.*t31.*t36.*t37+I2.*l.*r.*t2.*t4.*t15.*t17.*t31.*t36.*t38-I3.*l.*r.*t2.*t4.*t15.*t17.*t31.*t36.*t37-I3.*l.*r.*t2.*t4.*t15.*t17.*t31.*t36.*t38+I2.*l.*r.*t2.*t4.*t17.*t30.*t31.*u1.*u3-I3.*l.*r.*t2.*t4.*t17.*t30.*t31.*u1.*u3+I1.*l.*r.*t2.*t4.*t31.*t34.*t36.*u1.*u2.*2.0-I2.*l.*r.*t2.*t4.*t31.*t34.*t36.*u1.*u2.*3.0+I3.*l.*r.*t2.*t4.*t31.*t34.*t36.*u1.*u2+I1.*l.*r.*t2.*t4.*t17.*t30.*t31.*t34.*u1.*u3.*2.0+I2.*l.*r.*t2.*t4.*t15.*t30.*t31.*t34.*u2.*u3.*2.0-I3.*l.*r.*t2.*t4.*t15.*t30.*t31.*t34.*u2.*u3.*2.0-I3.*l.*r.*t2.*t4.*t17.*t30.*t31.*t34.*u1.*u3.*2.0)-t53.*(I2.*I3.*tau1p1+t4.*t22.*tau1p1+I3.*t6.*u2p1.*u3p1-I2.*t35.*u2p1.*u3p1+t2.*t4.*t5.*tau1p1+I3.*M.*t2.*tau1p1+I2.*M.*t5.*tau1p1+I3.*M.*t5.*tau1p1-I2.*I3.*f2p1.*r.*t3+I2.*M.*t2.*t11.*tau1p1-I3.*M.*t2.*t11.*tau1p1+I2.*t4.*t22.*u2p1.*u3p1-I3.*t4.*t22.*u2p1.*u3p1+M.*t5.*t6.*u2p1.*u3p1-M.*t2.*t35.*u2p1.*u3p1-M.*t5.*t35.*u2p1.*u3p1+t4.*t9.*t11.*t12.*tau1p1+l.*r.*t2.*t3.*t4.*tau3p1+l.*r.*t3.*t4.*t5.*tau3p1+t2.*t4.*t5.*t11.*t12.*tau1p1.*5.0+t3.*t4.*t7.*t8.*t9.*tau3p1-t4.*t8.*t9.*t10.*t11.*tau2p1+I2.*I3.*M.*t2.*u2p1.*u3p1-I2.*I3.*f3p1.*r.*t7.*t10-I3.*M.*f2p1.*r.*t2.*t3-I3.*M.*f2p1.*r.*t3.*t5+I2.*M.*l.*r.*t3.*tau3p1+I3.*M.*t2.*t11.*t12.*tau1p1+I2.*t2.*t4.*t5.*u2p1.*u3p1-I3.*t2.*t4.*t5.*u2p1.*u3p1+M.*t2.*t6.*t11.*u2p1.*u3p1+M.*t2.*t11.*t35.*u2p1.*u3p1-I2.*I3.*M.*t2.*t11.*u2p1.*u3p1.*2.0-I2.*M.*f2p1.*r.*t2.*t3.*t11+I3.*M.*f2p1.*r.*t2.*t3.*t11-I3.*M.*f3p1.*r.*t2.*t7.*t10-I2.*M.*f3p1.*r.*t5.*t7.*t10+I2.*M.*l.*r.*t7.*t8.*tau1p1.*2.0+I3.*M.*l.*r.*t7.*t8.*tau1p1.*2.0-I3.*M.*l.*r.*t7.*t10.*tau2p1+I2.*M.*t2.*t3.*t7.*t8.*tau3p1-I3.*M.*t2.*t8.*t10.*t11.*tau2p1-M.*l.*r.*t3.*t6.*u1p1.*u2p1+I2.*t4.*t9.*t11.*t12.*u2p1.*u3p1-I3.*t4.*t9.*t11.*t12.*u2p1.*u3p1-M.*t2.*t11.*t12.*t35.*u2p1.*u3p1+l.*r.*t2.*t4.*t7.*t8.*tau1p1.*2.0+l.*r.*t4.*t5.*t7.*t8.*tau1p1.*4.0-l.*r.*t2.*t4.*t7.*t10.*tau2p1-l.*r.*t4.*t5.*t7.*t10.*tau2p1+t2.*t3.*t4.*t5.*t7.*t8.*tau3p1.*3.0-t2.*t4.*t5.*t8.*t10.*t11.*tau2p1.*3.0+I1.*I2.*M.*l.*r.*t3.*u1p1.*u2p1+I2.*I3.*M.*l.*r.*t3.*u1p1.*u2p1+I2.*I3.*M.*t2.*t11.*t12.*u2p1.*u3p1+I2.*M.*f1p1.*l.*t2.*t3.*t7.*t10-I3.*M.*f1p1.*l.*t2.*t3.*t7.*t10-I3.*M.*f2p1.*l.*t2.*t3.*t7.*t8.*2.0-I2.*M.*f3p1.*l.*t2.*t8.*t10.*t11.*2.0+I2.*M.*f2p1.*r.*t2.*t3.*t11.*t12-I3.*M.*f2p1.*r.*t2.*t3.*t11.*t12-I2.*M.*f3p1.*r.*t2.*t7.*t10.*t11+I3.*M.*f3p1.*r.*t2.*t7.*t10.*t11+I2.*g.*r.*t3.*t4.*t5.*t7.*t10-I3.*g.*r.*t3.*t4.*t5.*t7.*t10+I1.*l.*r.*t2.*t3.*t4.*u1p1.*u2p1-I2.*l.*r.*t2.*t3.*t4.*u1p1.*u2p1+I3.*l.*r.*t2.*t3.*t4.*u1p1.*u2p1+I1.*l.*r.*t3.*t4.*t5.*u1p1.*u2p1-I2.*l.*r.*t3.*t4.*t5.*u1p1.*u2p1+I3.*l.*r.*t3.*t4.*t5.*u1p1.*u2p1+M.*l.*r.*t6.*t7.*t8.*u2p1.*u3p1.*2.0-M.*l.*r.*t7.*t10.*t35.*u1p1.*u3p1-M.*l.*r.*t7.*t8.*t35.*u2p1.*u3p1.*2.0+I2.*t2.*t3.*t4.*t5.*t7.*t10.*t54+I2.*t2.*t3.*t4.*t5.*t7.*t10.*t55-I3.*t2.*t3.*t4.*t5.*t7.*t10.*t54-I3.*t2.*t3.*t4.*t5.*t7.*t10.*t55+I1.*t3.*t4.*t7.*t8.*t9.*u1p1.*u2p1-I2.*t3.*t4.*t7.*t8.*t9.*u1p1.*u2p1+I1.*t4.*t8.*t9.*t10.*t11.*u1p1.*u3p1-I3.*t4.*t8.*t9.*t10.*t11.*u1p1.*u3p1+I2.*t2.*t4.*t5.*t11.*t12.*u2p1.*u3p1.*5.0-I3.*t2.*t4.*t5.*t11.*t12.*u2p1.*u3p1.*5.0-M.*t2.*t3.*t6.*t7.*t8.*u1p1.*u2p1-M.*t2.*t8.*t10.*t11.*t35.*u1p1.*u3p1+l.*r.*t2.*t3.*t4.*t11.*t12.*tau3p1.*2.0+I2.*M.*f1p1.*r.*t2.*t3.*t8.*t10.*t11-I3.*M.*f1p1.*r.*t2.*t3.*t8.*t10.*t11+I2.*g.*l.*t2.*t3.*t4.*t8.*t10.*t11-I3.*g.*l.*t2.*t3.*t4.*t8.*t10.*t11+I2.*l.*r.*t2.*t3.*t4.*t11.*u1p1.*u2p1-I3.*l.*r.*t2.*t3.*t4.*t11.*u1p1.*u2p1+I1.*l.*r.*t2.*t4.*t7.*t10.*u1p1.*u3p1+I1.*l.*r.*t4.*t5.*t7.*t10.*u1p1.*u3p1+I2.*l.*r.*t4.*t5.*t7.*t10.*u1p1.*u3p1-I3.*l.*r.*t4.*t5.*t7.*t10.*u1p1.*u3p1+I2.*l.*r.*t2.*t4.*t7.*t8.*u2p1.*u3p1.*2.0-I3.*l.*r.*t2.*t4.*t7.*t8.*u2p1.*u3p1.*2.0+I2.*l.*r.*t4.*t5.*t7.*t8.*u2p1.*u3p1.*4.0-I3.*l.*r.*t4.*t5.*t7.*t8.*u2p1.*u3p1.*4.0+I1.*t2.*t3.*t4.*t5.*t7.*t8.*u1p1.*u2p1.*3.0-I2.*t2.*t3.*t4.*t5.*t7.*t8.*u1p1.*u2p1.*3.0+I3.*t2.*t3.*t4.*t5.*t7.*t8.*u1p1.*u2p1.*2.0+I1.*t2.*t4.*t5.*t8.*t10.*t11.*u1p1.*u3p1.*3.0+I2.*t2.*t4.*t5.*t8.*t10.*t11.*u1p1.*u3p1.*2.0-I3.*t2.*t4.*t5.*t8.*t10.*t11.*u1p1.*u3p1.*3.0+l.*r.*t2.*t4.*t7.*t8.*t11.*t12.*tau1p1.*2.0-l.*r.*t2.*t4.*t7.*t10.*t11.*t12.*tau2p1.*2.0+I1.*I3.*M.*l.*r.*t7.*t10.*u1p1.*u3p1+I2.*I3.*M.*l.*r.*t7.*t10.*u1p1.*u3p1+I1.*I2.*M.*t2.*t3.*t7.*t8.*u1p1.*u2p1+I1.*I3.*M.*t2.*t8.*t10.*t11.*u1p1.*u3p1+I2.*l.*r.*t2.*t3.*t4.*t8.*t10.*t11.*t54+I2.*l.*r.*t2.*t3.*t4.*t8.*t10.*t11.*t55-I3.*l.*r.*t2.*t3.*t4.*t8.*t10.*t11.*t54-I3.*l.*r.*t2.*t3.*t4.*t8.*t10.*t11.*t55+I1.*l.*r.*t2.*t3.*t4.*t11.*t12.*u1p1.*u2p1.*2.0-I2.*l.*r.*t2.*t3.*t4.*t11.*t12.*u1p1.*u2p1.*3.0+I3.*l.*r.*t2.*t3.*t4.*t11.*t12.*u1p1.*u2p1+I2.*l.*r.*t2.*t4.*t7.*t10.*t11.*u1p1.*u3p1-I3.*l.*r.*t2.*t4.*t7.*t10.*t11.*u1p1.*u3p1+I1.*l.*r.*t2.*t4.*t7.*t10.*t11.*t12.*u1p1.*u3p1.*2.0-I3.*l.*r.*t2.*t4.*t7.*t10.*t11.*t12.*u1p1.*u3p1.*2.0+I2.*l.*r.*t2.*t4.*t7.*t8.*t11.*t12.*u2p1.*u3p1.*2.0-I3.*l.*r.*t2.*t4.*t7.*t8.*t11.*t12.*u2p1.*u3p1.*2.0)).*(1.0./8.0)).*(1.0./8.0)+dt.*t96.*(u2.*(1.0./2.0)+u2p1.*(1.0./2.0)+dt.*(t82.*(I1.*I3.*tau2+I1.*t35.*u1.*u3-I3.*t83.*u1.*u3+t4.*t9.*t31.*tau2-I1.*I3.*f3.*l+I3.*M.*t2.*tau2+I1.*M.*t5.*tau2-I3.*M.*f3.*l.*t2-I1.*M.*f3.*l.*t5+I1.*I3.*f1.*r.*t36+I1.*M.*t2.*t31.*tau2+I1.*t4.*t22.*u1.*u3+M.*t2.*t35.*u1.*u3-M.*t5.*t83.*u1.*u3+t2.*t4.*t5.*t31.*tau2-t4.*t9.*t31.*t34.*tau2-t4.*t9.*t15.*t17.*t31.*tau1-t2.*t4.*t5.*t31.*t34.*tau2-t4.*t9.*t17.*t30.*t36.*tau3+I1.*I3.*M.*g.*l.*t36-I1.*I3.*M.*t2.*u1.*u3+I1.*I3.*M.*t5.*u1.*u3.*2.0-I1.*M.*f3.*l.*t2.*t31+I3.*M.*f3.*l.*t2.*t31-I1.*I3.*f3.*r.*t15.*t30+I3.*M.*f1.*r.*t2.*t36+I1.*M.*f1.*r.*t5.*t36-I3.*M.*t2.*t31.*t34.*tau2+I3.*g.*l.*t2.*t4.*t36+I1.*g.*l.*t4.*t5.*t36+I3.*t2.*t4.*t5.*u1.*u3-I1.*t4.*t9.*t31.*u1.*u3+I3.*t4.*t9.*t31.*u1.*u3-M.*t2.*t31.*t83.*u1.*u3+I1.*I3.*M.*l.*r.*t36.*t37+I1.*I3.*M.*l.*r.*t36.*t38+I1.*I3.*M.*t2.*t31.*u1.*u3-I1.*M.*f3.*l.*t2.*t31.*t34.*2.0-I3.*M.*f3.*r.*t2.*t15.*t30-I1.*M.*f3.*r.*t5.*t15.*t30.*3.0+I1.*M.*l.*r.*t15.*t30.*tau2.*2.0-I3.*M.*l.*r.*t17.*t30.*tau1-I3.*M.*t2.*t15.*t17.*t31.*tau1-I1.*M.*t2.*t17.*t30.*t36.*tau3+I3.*l.*r.*t2.*t4.*t36.*t37+I1.*l.*r.*t4.*t5.*t36.*t37+I3.*l.*r.*t2.*t4.*t36.*t38+I1.*l.*r.*t4.*t5.*t36.*t38+I1.*t4.*t9.*t31.*t34.*u1.*u3-I3.*t4.*t9.*t31.*t34.*u1.*u3-M.*t2.*t31.*t34.*t35.*u1.*u3-l.*r.*t2.*t4.*t17.*t30.*tau1-l.*r.*t4.*t5.*t17.*t30.*tau1-t2.*t4.*t5.*t15.*t17.*t31.*tau1.*3.0-t2.*t4.*t5.*t17.*t30.*t36.*tau3+I1.*I3.*M.*t2.*t31.*t34.*u1.*u3+I1.*M.*f1.*l.*t2.*t15.*t30.*t36.*2.0-I1.*M.*f2.*l.*t2.*t17.*t30.*t36+I3.*M.*f2.*l.*t2.*t17.*t30.*t36-I1.*M.*f3.*r.*t2.*t15.*t30.*t31+I3.*M.*f3.*r.*t2.*t15.*t30.*t31+I1.*M.*f1.*r.*t2.*t31.*t34.*t36-I3.*M.*f1.*r.*t2.*t31.*t34.*t36+I1.*g.*l.*t2.*t4.*t31.*t34.*t36-I3.*g.*l.*t2.*t4.*t31.*t34.*t36+I1.*g.*r.*t4.*t5.*t15.*t30.*t36.*2.0+M.*l.*r.*t17.*t30.*t35.*u2.*u3-M.*l.*r.*t15.*t30.*t83.*u1.*u3.*2.0+I1.*t2.*t4.*t5.*t15.*t30.*t36.*t37.*2.0+I1.*t2.*t4.*t5.*t15.*t30.*t36.*t38.*2.0+I1.*t2.*t4.*t5.*t31.*t34.*u1.*u3.*3.0-I2.*t4.*t9.*t15.*t17.*t31.*u2.*u3-I3.*t2.*t4.*t5.*t31.*t34.*u1.*u3+I3.*t4.*t9.*t15.*t17.*t31.*u2.*u3-I1.*t4.*t9.*t17.*t30.*t36.*u1.*u2+I2.*t4.*t9.*t17.*t30.*t36.*u1.*u2+M.*t2.*t15.*t17.*t31.*t35.*u2.*u3-M.*t2.*t17.*t30.*t36.*t83.*u1.*u2+l.*r.*t2.*t4.*t15.*t30.*t31.*tau2.*2.0-I1.*M.*f2.*r.*t2.*t15.*t17.*t31.*t36+I3.*M.*f2.*r.*t2.*t15.*t17.*t31.*t36+I1.*l.*r.*t2.*t4.*t31.*t34.*t36.*t37+I1.*l.*r.*t2.*t4.*t31.*t34.*t36.*t38-I3.*l.*r.*t2.*t4.*t31.*t34.*t36.*t37-I3.*l.*r.*t2.*t4.*t31.*t34.*t36.*t38+I3.*l.*r.*t2.*t4.*t15.*t30.*u1.*u3+I1.*l.*r.*t4.*t5.*t15.*t30.*u1.*u3.*3.0-I2.*l.*r.*t2.*t4.*t17.*t30.*u2.*u3+I3.*l.*r.*t2.*t4.*t17.*t30.*u2.*u3-I2.*l.*r.*t4.*t5.*t17.*t30.*u2.*u3+I3.*l.*r.*t4.*t5.*t17.*t30.*u2.*u3-I2.*t2.*t4.*t5.*t15.*t17.*t31.*u2.*u3.*3.0+I3.*t2.*t4.*t5.*t15.*t17.*t31.*u2.*u3.*3.0+I2.*t2.*t4.*t5.*t17.*t30.*t36.*u1.*u2-I3.*t2.*t4.*t5.*t17.*t30.*t36.*u1.*u2-l.*r.*t2.*t4.*t15.*t17.*t31.*t36.*tau3.*2.0-l.*r.*t2.*t4.*t15.*t30.*t31.*t34.*tau2.*2.0-l.*r.*t2.*t4.*t17.*t30.*t31.*t34.*tau1.*2.0+I1.*I3.*M.*l.*r.*t15.*t30.*u1.*u3.*3.0-I2.*I3.*M.*l.*r.*t17.*t30.*u2.*u3-I2.*I3.*M.*t2.*t15.*t17.*t31.*u2.*u3+I1.*I2.*M.*t2.*t17.*t30.*t36.*u1.*u2-I1.*l.*r.*t2.*t4.*t15.*t30.*t31.*u1.*u3+I3.*l.*r.*t2.*t4.*t15.*t30.*t31.*u1.*u3-I1.*l.*r.*t2.*t4.*t15.*t17.*t31.*t36.*u1.*u2+I2.*l.*r.*t2.*t4.*t15.*t17.*t31.*t36.*u1.*u2.*2.0-I3.*l.*r.*t2.*t4.*t15.*t17.*t31.*t36.*u1.*u2+I1.*l.*r.*t2.*t4.*t15.*t30.*t31.*t34.*u1.*u3.*2.0-I3.*l.*r.*t2.*t4.*t15.*t30.*t31.*t34.*u1.*u3.*2.0-I2.*l.*r.*t2.*t4.*t17.*t30.*t31.*t34.*u2.*u3.*2.0+I3.*l.*r.*t2.*t4.*t17.*t30.*t31.*t34.*u2.*u3.*2.0)-t53.*(I1.*I3.*tau2p1+I1.*t35.*u1p1.*u3p1-I3.*t83.*u1p1.*u3p1+t4.*t9.*t11.*tau2p1-I1.*I3.*f3p1.*l+I3.*M.*t2.*tau2p1+I1.*M.*t5.*tau2p1-I3.*M.*f3p1.*l.*t2-I1.*M.*f3p1.*l.*t5+I1.*I3.*f1p1.*r.*t3+I1.*M.*t2.*t11.*tau2p1+I1.*t4.*t22.*u1p1.*u3p1+M.*t2.*t35.*u1p1.*u3p1-M.*t5.*t83.*u1p1.*u3p1+t2.*t4.*t5.*t11.*tau2p1-t4.*t9.*t11.*t12.*tau2p1-t4.*t8.*t9.*t10.*t11.*tau1p1-t2.*t4.*t5.*t11.*t12.*tau2p1-t3.*t4.*t7.*t9.*t10.*tau3p1+I1.*I3.*M.*g.*l.*t3-I1.*I3.*M.*t2.*u1p1.*u3p1+I1.*I3.*M.*t5.*u1p1.*u3p1.*2.0-I1.*M.*f3p1.*l.*t2.*t11+I3.*M.*f3p1.*l.*t2.*t11-I1.*I3.*f3p1.*r.*t7.*t8+I3.*M.*f1p1.*r.*t2.*t3+I1.*M.*f1p1.*r.*t3.*t5-I3.*M.*t2.*t11.*t12.*tau2p1+I3.*g.*l.*t2.*t3.*t4+I1.*g.*l.*t3.*t4.*t5+I3.*t2.*t4.*t5.*u1p1.*u3p1-I1.*t4.*t9.*t11.*u1p1.*u3p1+I3.*t4.*t9.*t11.*u1p1.*u3p1-M.*t2.*t11.*t83.*u1p1.*u3p1+I1.*I3.*M.*l.*r.*t3.*t54+I1.*I3.*M.*l.*r.*t3.*t55+I1.*I3.*M.*t2.*t11.*u1p1.*u3p1-I1.*M.*f3p1.*l.*t2.*t11.*t12.*2.0-I3.*M.*f3p1.*r.*t2.*t7.*t8-I1.*M.*f3p1.*r.*t5.*t7.*t8.*3.0-I3.*M.*l.*r.*t7.*t10.*tau1p1+I1.*M.*l.*r.*t7.*t8.*tau2p1.*2.0-I3.*M.*t2.*t8.*t10.*t11.*tau1p1-I1.*M.*t2.*t3.*t7.*t10.*tau3p1+I3.*l.*r.*t2.*t3.*t4.*t54+I1.*l.*r.*t3.*t4.*t5.*t54+I3.*l.*r.*t2.*t3.*t4.*t55+I1.*l.*r.*t3.*t4.*t5.*t55+I1.*t4.*t9.*t11.*t12.*u1p1.*u3p1-I3.*t4.*t9.*t11.*t12.*u1p1.*u3p1-M.*t2.*t11.*t12.*t35.*u1p1.*u3p1-l.*r.*t2.*t4.*t7.*t10.*tau1p1-l.*r.*t4.*t5.*t7.*t10.*tau1p1-t2.*t4.*t5.*t8.*t10.*t11.*tau1p1.*3.0-t2.*t3.*t4.*t5.*t7.*t10.*tau3p1+I1.*I3.*M.*t2.*t11.*t12.*u1p1.*u3p1+I1.*M.*f1p1.*l.*t2.*t3.*t7.*t8.*2.0-I1.*M.*f2p1.*l.*t2.*t3.*t7.*t10+I3.*M.*f2p1.*l.*t2.*t3.*t7.*t10+I1.*M.*f1p1.*r.*t2.*t3.*t11.*t12-I3.*M.*f1p1.*r.*t2.*t3.*t11.*t12-I1.*M.*f3p1.*r.*t2.*t7.*t8.*t11+I3.*M.*f3p1.*r.*t2.*t7.*t8.*t11+I1.*g.*l.*t2.*t3.*t4.*t11.*t12-I3.*g.*l.*t2.*t3.*t4.*t11.*t12+I1.*g.*r.*t3.*t4.*t5.*t7.*t8.*2.0+M.*l.*r.*t7.*t10.*t35.*u2p1.*u3p1-M.*l.*r.*t7.*t8.*t83.*u1p1.*u3p1.*2.0+I1.*t2.*t3.*t4.*t5.*t7.*t8.*t54.*2.0+I1.*t2.*t3.*t4.*t5.*t7.*t8.*t55.*2.0-I1.*t3.*t4.*t7.*t9.*t10.*u1p1.*u2p1+I2.*t3.*t4.*t7.*t9.*t10.*u1p1.*u2p1+I1.*t2.*t4.*t5.*t11.*t12.*u1p1.*u3p1.*3.0-I3.*t2.*t4.*t5.*t11.*t12.*u1p1.*u3p1-I2.*t4.*t8.*t9.*t10.*t11.*u2p1.*u3p1+I3.*t4.*t8.*t9.*t10.*t11.*u2p1.*u3p1+M.*t2.*t8.*t10.*t11.*t35.*u2p1.*u3p1-M.*t2.*t3.*t7.*t10.*t83.*u1p1.*u2p1+l.*r.*t2.*t4.*t7.*t8.*t11.*tau2p1.*2.0-I1.*M.*f2p1.*r.*t2.*t3.*t8.*t10.*t11+I3.*M.*f2p1.*r.*t2.*t3.*t8.*t10.*t11+I1.*l.*r.*t2.*t3.*t4.*t11.*t12.*t54+I1.*l.*r.*t2.*t3.*t4.*t11.*t12.*t55-I3.*l.*r.*t2.*t3.*t4.*t11.*t12.*t54-I3.*l.*r.*t2.*t3.*t4.*t11.*t12.*t55+I3.*l.*r.*t2.*t4.*t7.*t8.*u1p1.*u3p1+I1.*l.*r.*t4.*t5.*t7.*t8.*u1p1.*u3p1.*3.0-I2.*l.*r.*t2.*t4.*t7.*t10.*u2p1.*u3p1+I3.*l.*r.*t2.*t4.*t7.*t10.*u2p1.*u3p1-I2.*l.*r.*t4.*t5.*t7.*t10.*u2p1.*u3p1+I3.*l.*r.*t4.*t5.*t7.*t10.*u2p1.*u3p1+I2.*t2.*t3.*t4.*t5.*t7.*t10.*u1p1.*u2p1-I3.*t2.*t3.*t4.*t5.*t7.*t10.*u1p1.*u2p1-I2.*t2.*t4.*t5.*t8.*t10.*t11.*u2p1.*u3p1.*3.0+I3.*t2.*t4.*t5.*t8.*t10.*t11.*u2p1.*u3p1.*3.0-l.*r.*t2.*t4.*t7.*t10.*t11.*t12.*tau1p1.*2.0-l.*r.*t2.*t4.*t7.*t8.*t11.*t12.*tau2p1.*2.0-l.*r.*t2.*t3.*t4.*t8.*t10.*t11.*tau3p1.*2.0+I1.*I3.*M.*l.*r.*t7.*t8.*u1p1.*u3p1.*3.0-I2.*I3.*M.*l.*r.*t7.*t10.*u2p1.*u3p1+I1.*I2.*M.*t2.*t3.*t7.*t10.*u1p1.*u2p1-I2.*I3.*M.*t2.*t8.*t10.*t11.*u2p1.*u3p1-I1.*l.*r.*t2.*t4.*t7.*t8.*t11.*u1p1.*u3p1+I3.*l.*r.*t2.*t4.*t7.*t8.*t11.*u1p1.*u3p1-I1.*l.*r.*t2.*t3.*t4.*t8.*t10.*t11.*u1p1.*u2p1+I2.*l.*r.*t2.*t3.*t4.*t8.*t10.*t11.*u1p1.*u2p1.*2.0-I3.*l.*r.*t2.*t3.*t4.*t8.*t10.*t11.*u1p1.*u2p1+I1.*l.*r.*t2.*t4.*t7.*t8.*t11.*t12.*u1p1.*u3p1.*2.0-I3.*l.*r.*t2.*t4.*t7.*t8.*t11.*t12.*u1p1.*u3p1.*2.0-I2.*l.*r.*t2.*t4.*t7.*t10.*t11.*t12.*u2p1.*u3p1.*2.0+I3.*l.*r.*t2.*t4.*t7.*t10.*t11.*t12.*u2p1.*u3p1.*2.0)).*(1.0./8.0)).*(1.0./8.0)-dt.*t53.*t96.*(I3.*t6.*u2p1-I2.*t35.*u2p1+I2.*t4.*t22.*u2p1-I3.*t4.*t22.*u2p1+M.*t5.*t6.*u2p1-M.*t2.*t35.*u2p1-M.*t5.*t35.*u2p1+I2.*I3.*M.*t2.*u2p1+I2.*t2.*t4.*t5.*u2p1-I3.*t2.*t4.*t5.*u2p1+M.*t2.*t6.*t11.*u2p1+M.*t2.*t11.*t35.*u2p1-I2.*I3.*M.*t2.*t11.*u2p1.*2.0+I2.*t4.*t9.*t11.*t12.*u2p1-I3.*t4.*t9.*t11.*t12.*u2p1-M.*t2.*t11.*t12.*t35.*u2p1+I2.*I3.*M.*t2.*t11.*t12.*u2p1+M.*l.*r.*t6.*t7.*t8.*u2p1.*2.0-M.*l.*r.*t7.*t10.*t35.*u1p1-M.*l.*r.*t7.*t8.*t35.*u2p1.*2.0+I1.*t4.*t8.*t9.*t10.*t11.*u1p1-I3.*t4.*t8.*t9.*t10.*t11.*u1p1+I2.*t2.*t4.*t5.*t11.*t12.*u2p1.*5.0-I3.*t2.*t4.*t5.*t11.*t12.*u2p1.*5.0-M.*t2.*t8.*t10.*t11.*t35.*u1p1+I1.*I3.*M.*l.*r.*t7.*t10.*u1p1+I2.*I3.*M.*l.*r.*t7.*t10.*u1p1+I1.*I3.*M.*t2.*t8.*t10.*t11.*u1p1+I1.*l.*r.*t2.*t4.*t7.*t10.*u1p1+I1.*l.*r.*t4.*t5.*t7.*t10.*u1p1+I2.*l.*r.*t4.*t5.*t7.*t10.*u1p1-I3.*l.*r.*t4.*t5.*t7.*t10.*u1p1+I2.*l.*r.*t2.*t4.*t7.*t8.*u2p1.*2.0-I3.*l.*r.*t2.*t4.*t7.*t8.*u2p1.*2.0+I2.*l.*r.*t4.*t5.*t7.*t8.*u2p1.*4.0-I3.*l.*r.*t4.*t5.*t7.*t8.*u2p1.*4.0+I1.*t2.*t4.*t5.*t8.*t10.*t11.*u1p1.*3.0+I2.*t2.*t4.*t5.*t8.*t10.*t11.*u1p1.*2.0-I3.*t2.*t4.*t5.*t8.*t10.*t11.*u1p1.*3.0+I2.*t2.*t3.*t4.*t5.*t7.*t10.*u3p1.*2.0-I3.*t2.*t3.*t4.*t5.*t7.*t10.*u3p1.*2.0+I2.*l.*r.*t2.*t4.*t7.*t10.*t11.*u1p1-I3.*l.*r.*t2.*t4.*t7.*t10.*t11.*u1p1+I1.*l.*r.*t2.*t4.*t7.*t10.*t11.*t12.*u1p1.*2.0-I3.*l.*r.*t2.*t4.*t7.*t10.*t11.*t12.*u1p1.*2.0+I2.*l.*r.*t2.*t4.*t7.*t8.*t11.*t12.*u2p1.*2.0-I3.*l.*r.*t2.*t4.*t7.*t8.*t11.*t12.*u2p1.*2.0+I2.*l.*r.*t2.*t3.*t4.*t8.*t10.*t11.*u3p1.*2.0-I3.*l.*r.*t2.*t3.*t4.*t8.*t10.*t11.*u3p1.*2.0).*(1.0./8.0)+dt.*t53.*t95.*(I1.*t35.*u1p1-I3.*t83.*u1p1+I1.*t4.*t22.*u1p1+M.*t2.*t35.*u1p1-M.*t5.*t83.*u1p1-I1.*I3.*M.*t2.*u1p1+I1.*I3.*M.*t5.*u1p1.*2.0+I3.*t2.*t4.*t5.*u1p1-I1.*t4.*t9.*t11.*u1p1+I3.*t4.*t9.*t11.*u1p1-M.*t2.*t11.*t83.*u1p1+I1.*I3.*M.*t2.*t11.*u1p1+I1.*t4.*t9.*t11.*t12.*u1p1-I3.*t4.*t9.*t11.*t12.*u1p1-M.*t2.*t11.*t12.*t35.*u1p1+I1.*I3.*M.*l.*r.*t3.*u3p1.*2.0+I1.*I3.*M.*t2.*t11.*t12.*u1p1+I3.*l.*r.*t2.*t3.*t4.*u3p1.*2.0+I1.*l.*r.*t3.*t4.*t5.*u3p1.*2.0+M.*l.*r.*t7.*t10.*t35.*u2p1-M.*l.*r.*t7.*t8.*t83.*u1p1.*2.0+I1.*t2.*t4.*t5.*t11.*t12.*u1p1.*3.0-I3.*t2.*t4.*t5.*t11.*t12.*u1p1-I2.*t4.*t8.*t9.*t10.*t11.*u2p1+I3.*t4.*t8.*t9.*t10.*t11.*u2p1+M.*t2.*t8.*t10.*t11.*t35.*u2p1+I1.*I3.*M.*l.*r.*t7.*t8.*u1p1.*3.0-I2.*I3.*M.*l.*r.*t7.*t10.*u2p1-I2.*I3.*M.*t2.*t8.*t10.*t11.*u2p1+I3.*l.*r.*t2.*t4.*t7.*t8.*u1p1+I1.*l.*r.*t4.*t5.*t7.*t8.*u1p1.*3.0-I2.*l.*r.*t2.*t4.*t7.*t10.*u2p1+I3.*l.*r.*t2.*t4.*t7.*t10.*u2p1-I2.*l.*r.*t4.*t5.*t7.*t10.*u2p1+I3.*l.*r.*t4.*t5.*t7.*t10.*u2p1+I1.*t2.*t3.*t4.*t5.*t7.*t8.*u3p1.*4.0-I2.*t2.*t4.*t5.*t8.*t10.*t11.*u2p1.*3.0+I3.*t2.*t4.*t5.*t8.*t10.*t11.*u2p1.*3.0-I1.*l.*r.*t2.*t4.*t7.*t8.*t11.*u1p1+I3.*l.*r.*t2.*t4.*t7.*t8.*t11.*u1p1+I1.*l.*r.*t2.*t3.*t4.*t11.*t12.*u3p1.*2.0-I3.*l.*r.*t2.*t3.*t4.*t11.*t12.*u3p1.*2.0+I1.*l.*r.*t2.*t4.*t7.*t8.*t11.*t12.*u1p1.*2.0-I3.*l.*r.*t2.*t4.*t7.*t8.*t11.*t12.*u1p1.*2.0-I2.*l.*r.*t2.*t4.*t7.*t10.*t11.*t12.*u2p1.*2.0+I3.*l.*r.*t2.*t4.*t7.*t10.*t11.*t12.*u2p1.*2.0).*(1.0./8.0)))./cos(t21)-3.0./4.0;
