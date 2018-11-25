function out1 = dhi_dx_1_7(q1,q2,q3,q4,q5,u1,u2,u3,tau1,tau2,tau3,f1,f2,f3,q1p1,q2p1,q3p1,q4p1,q5p1,u1p1,u2p1,u3p1,tau1p1,tau2p1,tau3p1,f1p1,f2p1,f3p1,dt,M,l,r,g,I1,I2,I3)
%DHI_DX_1_7
%    OUT1 = DHI_DX_1_7(Q1,Q2,Q3,Q4,Q5,U1,U2,U3,TAU1,TAU2,TAU3,F1,F2,F3,Q1P1,Q2P1,Q3P1,Q4P1,Q5P1,U1P1,U2P1,U3P1,TAU1P1,TAU2P1,TAU3P1,F1P1,F2P1,F3P1,DT,M,L,R,G,I1,I2,I3)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    23-Nov-2018 20:37:19

t2 = cos(q3);
t3 = sin(q3);
t4 = cos(q3p1);
t5 = sin(q3p1);
t6 = l.^2;
t7 = sin(q2);
t8 = M.^2;
t9 = r.^2;
t10 = cos(q2);
t11 = t9.^2;
t12 = I3.^2;
t13 = t10.^2;
t14 = t2.^2;
t15 = q3.*(1.0./2.0);
t16 = q3p1.*(1.0./2.0);
t17 = 1.0./t10;
t18 = t2.*u1;
t55 = t3.*u2;
t19 = t18-t55;
t20 = cos(q2p1);
t21 = 1.0./t20;
t22 = sin(q2p1);
t23 = t4.*u1p1;
t57 = t5.*u2p1;
t24 = t23-t57;
t25 = t21.*t22.*t24;
t56 = t7.*t17.*t19;
t26 = t25-t56+u3-u3p1;
t27 = dt.*t26.*(1.0./8.0);
t28 = t15+t16+t27;
t29 = t6.^2;
t30 = I1.*t8.*t29;
t31 = I3.*t8.*t11;
t32 = I1.*I2.*I3;
t33 = I1.*I2.*M.*t6;
t34 = I1.*I3.*M.*t6;
t35 = I1.*I3.*M.*t9;
t36 = I2.*I3.*M.*t9;
t37 = I2.*t8.*t11.*t13;
t38 = I1.*t6.*t8.*t9;
t39 = I3.*t6.*t8.*t9;
t40 = I2.*t6.*t8.*t9.*t13;
t41 = I1.*I2.*M.*t9.*t13;
t42 = I1.*t8.*t11.*t13.*t14;
t43 = I1.*l.*r.*t2.*t8.*t9.*t10.*2.0;
t44 = I1.*l.*r.*t2.*t6.*t8.*t10.*4.0;
t45 = I3.*l.*r.*t2.*t8.*t9.*t10.*2.0;
t46 = I1.*t6.*t8.*t9.*t13.*t14.*5.0;
t47 = I1.*I3.*M.*t9.*t13.*t14;
t48 = I2.*l.*r.*t2.*t8.*t9.*t10.*t13.*2.0;
t49 = I1.*l.*r.*t2.*t8.*t9.*t10.*t13.*t14.*2.0;
t50 = I1.*I2.*M.*l.*r.*t2.*t10.*2.0;
t51 = I1.*I3.*M.*l.*r.*t2.*t10.*2.0;
t59 = I3.*t8.*t11.*t13;
t60 = I3.*t6.*t8.*t9.*t13;
t61 = I1.*I3.*M.*t9.*t13;
t62 = I2.*t8.*t11.*t13.*t14;
t63 = I2.*t6.*t8.*t9.*t13.*t14;
t64 = I2.*I3.*M.*t9.*t13.*t14;
t65 = I3.*l.*r.*t2.*t8.*t9.*t10.*t13.*2.0;
t66 = I2.*l.*r.*t2.*t8.*t9.*t10.*t13.*t14.*2.0;
t52 = t30+t31+t32+t33+t34+t35+t36+t37+t38+t39+t40+t41+t42+t43+t44+t45+t46+t47+t48+t49+t50+t51-t59-t60-t61-t62-t63-t64-t65-t66;
t53 = 1.0./t52;
t54 = I2.^2;
t58 = cos(t28);
t67 = I1.^2;
t68 = u2.^2;
t69 = u3.^2;
t70 = t20.^2;
t71 = t4.^2;
t72 = u2p1.^2;
t73 = u3p1.^2;
t74 = sin(t28);
t75 = I2.*t8.*t11.*t70;
t76 = I2.*t6.*t8.*t9.*t70;
t77 = I1.*I2.*M.*t9.*t70;
t78 = I1.*t8.*t11.*t70.*t71;
t79 = I1.*l.*r.*t4.*t8.*t9.*t20.*2.0;
t80 = I1.*l.*r.*t4.*t6.*t8.*t20.*4.0;
t81 = I3.*l.*r.*t4.*t8.*t9.*t20.*2.0;
t82 = I1.*t6.*t8.*t9.*t70.*t71.*5.0;
t83 = I1.*I3.*M.*t9.*t70.*t71;
t84 = I2.*l.*r.*t4.*t8.*t9.*t20.*t70.*2.0;
t85 = I1.*l.*r.*t4.*t8.*t9.*t20.*t70.*t71.*2.0;
t86 = I1.*I2.*M.*l.*r.*t4.*t20.*2.0;
t87 = I1.*I3.*M.*l.*r.*t4.*t20.*2.0;
t161 = I3.*t8.*t11.*t70;
t162 = I3.*t6.*t8.*t9.*t70;
t163 = I1.*I3.*M.*t9.*t70;
t164 = I2.*t8.*t11.*t70.*t71;
t165 = I2.*t6.*t8.*t9.*t70.*t71;
t166 = I2.*I3.*M.*t9.*t70.*t71;
t167 = I3.*l.*r.*t4.*t8.*t9.*t20.*t70.*2.0;
t168 = I2.*l.*r.*t4.*t8.*t9.*t20.*t70.*t71.*2.0;
t88 = t30+t31+t32+t33+t34+t35+t36+t38+t39+t75+t76+t77+t78+t79+t80+t81+t82+t83+t84+t85+t86+t87-t161-t162-t163-t164-t165-t166-t167-t168;
t89 = 1.0./t88;
t90 = q2.*(1.0./2.0);
t91 = q2p1.*(1.0./2.0);
t92 = t2.*u2;
t93 = t3.*u1;
t98 = t4.*u2p1;
t99 = t5.*u1p1;
t94 = t92+t93-t98-t99;
t95 = dt.*t94.*(1.0./8.0);
t96 = t90+t91+t95;
t97 = cos(t96);
t100 = u2.*(1.0./2.0);
t101 = u2p1.*(1.0./2.0);
t102 = I1.*I3.*tau2;
t103 = I1.*M.*t6.*tau2;
t104 = I1.*t12.*u1.*u3;
t105 = I3.*M.*t9.*tau2;
t106 = t8.*t11.*t13.*tau2;
t107 = I1.*t8.*t29.*u1.*u3;
t108 = M.*t9.*t12.*u1.*u3;
t109 = I3.*M.*f1.*r.*t7.*t9;
t110 = t6.*t8.*t9.*t13.*tau2;
t111 = I1.*g.*l.*t6.*t7.*t8;
t112 = I1.*M.*t9.*t13.*tau2;
t113 = I1.*I3.*f1.*r.*t7;
t114 = I3.*M.*f3.*r.*t2.*t9.*t10.*t13;
t115 = I3.*t8.*t11.*t13.*u1.*u3;
t116 = I1.*l.*r.*t6.*t7.*t8.*t68;
t117 = I1.*l.*r.*t6.*t7.*t8.*t69;
t118 = I3.*l.*r.*t7.*t8.*t9.*t68;
t119 = I3.*l.*r.*t7.*t8.*t9.*t69;
t120 = I1.*M.*f1.*r.*t6.*t7;
t121 = l.*r.*t2.*t8.*t9.*t10.*t13.*tau2.*2.0;
t122 = I3.*M.*f3.*l.*t9.*t13;
t123 = I3.*t6.*t8.*t9.*u1.*u3;
t124 = I3.*g.*l.*t7.*t8.*t9;
t125 = I1.*I3.*M.*t6.*u1.*u3.*2.0;
t126 = I1.*I3.*M.*g.*l.*t7;
t127 = I1.*I3.*M.*t9.*t13.*u1.*u3;
t128 = I1.*t8.*t11.*t13.*t14.*u1.*u3;
t129 = I1.*M.*f1.*r.*t7.*t9.*t13.*t14;
t130 = I1.*I3.*M.*l.*r.*t7.*t68;
t131 = I1.*I3.*M.*l.*r.*t7.*t69;
t132 = I1.*M.*l.*r.*t2.*t10.*tau2.*2.0;
t133 = I1.*l.*r.*t2.*t6.*t8.*t10.*u1.*u3.*3.0;
t134 = I3.*l.*r.*t2.*t8.*t9.*t10.*u1.*u3;
t135 = I1.*M.*f1.*l.*t2.*t7.*t9.*t10.*2.0;
t136 = I3.*l.*r.*t3.*t8.*t9.*t10.*u2.*u3;
t137 = I3.*l.*r.*t3.*t6.*t8.*t10.*u2.*u3;
t138 = I3.*M.*f2.*l.*t3.*t7.*t9.*t10;
t139 = I1.*t6.*t8.*t9.*t13.*t14.*u1.*u3.*3.0;
t140 = I1.*g.*l.*t7.*t8.*t9.*t13.*t14;
t141 = I1.*t2.*t6.*t7.*t8.*t9.*t10.*t68.*2.0;
t142 = I1.*t2.*t6.*t7.*t8.*t9.*t10.*t69.*2.0;
t143 = I1.*I3.*M.*t9.*t13.*t14.*u1.*u3;
t144 = I3.*l.*r.*t2.*t8.*t9.*t10.*t13.*u1.*u3;
t145 = I1.*g.*r.*t2.*t6.*t7.*t8.*t10.*2.0;
t146 = I2.*t3.*t7.*t8.*t10.*t11.*u1.*u2;
t147 = I1.*l.*r.*t7.*t8.*t9.*t13.*t14.*t68;
t148 = I1.*l.*r.*t7.*t8.*t9.*t13.*t14.*t69;
t149 = M.*l.*r.*t3.*t10.*t12.*u2.*u3;
t150 = I1.*l.*r.*t2.*t8.*t9.*t10.*t13.*t14.*u1.*u3.*2.0;
t151 = M.*t2.*t3.*t9.*t12.*t13.*u2.*u3;
t152 = I3.*t2.*t3.*t8.*t11.*t13.*u2.*u3;
t153 = I3.*M.*f2.*r.*t2.*t3.*t7.*t9.*t13;
t154 = I2.*t3.*t6.*t7.*t8.*t9.*t10.*u1.*u2;
t155 = I1.*I2.*M.*t3.*t7.*t9.*t10.*u1.*u2;
t156 = I1.*I3.*M.*l.*r.*t2.*t10.*u1.*u3.*3.0;
t157 = I3.*t2.*t3.*t6.*t8.*t9.*t13.*u2.*u3.*3.0;
t158 = I3.*l.*r.*t3.*t8.*t9.*t10.*t13.*t14.*u2.*u3.*2.0;
t159 = I2.*l.*r.*t2.*t3.*t7.*t8.*t9.*t13.*u1.*u2.*2.0;
t160 = t102+t103+t104+t105+t106+t107+t108+t109+t110+t111+t112+t113+t114+t115+t116+t117+t118+t119+t120+t121+t122+t123+t124+t125+t126+t127+t128+t129+t130+t131+t132+t133+t134+t135+t136+t137+t138+t139+t140+t141+t142+t143+t144+t145+t146+t147+t148+t149+t150+t151+t152+t153+t154+t155+t156+t157+t158+t159-I3.*t67.*u1.*u3-I1.*I3.*f3.*l-I1.*M.*f3.*l.*t6-I3.*M.*f3.*l.*t9-M.*t6.*t67.*u1.*u3-t8.*t11.*t13.*t14.*tau2-t2.*t3.*t8.*t11.*t13.*tau1-t3.*t7.*t8.*t10.*t11.*tau3-t6.*t8.*t9.*t13.*t14.*tau2-I1.*I3.*M.*t9.*u1.*u3-I1.*M.*f3.*l.*t9.*t13-I1.*I3.*f3.*r.*t2.*t10-I3.*M.*t9.*t13.*t14.*tau2-I1.*t8.*t11.*t13.*u1.*u3-M.*t9.*t13.*t67.*u1.*u3-I1.*M.*f3.*l.*t9.*t13.*t14.*2.0-I1.*M.*f3.*r.*t2.*t6.*t10.*3.0-I3.*M.*f3.*r.*t2.*t9.*t10-I3.*M.*l.*r.*t3.*t10.*tau1-I3.*M.*t2.*t3.*t9.*t13.*tau1-I1.*M.*t3.*t7.*t9.*t10.*tau3-I3.*t8.*t11.*t13.*t14.*u1.*u3-M.*t9.*t12.*t13.*t14.*u1.*u3-l.*r.*t3.*t6.*t8.*t10.*tau1-l.*r.*t3.*t8.*t9.*t10.*tau1-t2.*t3.*t6.*t8.*t9.*t13.*tau1.*3.0-t3.*t6.*t7.*t8.*t9.*t10.*tau3-I1.*M.*f2.*l.*t3.*t7.*t9.*t10-I1.*M.*f3.*r.*t2.*t9.*t10.*t13-I3.*M.*f1.*r.*t7.*t9.*t13.*t14-I3.*g.*l.*t7.*t8.*t9.*t13.*t14-M.*l.*r.*t2.*t10.*t67.*u1.*u3.*2.0-I1.*t3.*t7.*t8.*t10.*t11.*u1.*u2-I2.*t2.*t3.*t8.*t11.*t13.*u2.*u3-I3.*t6.*t8.*t9.*t13.*t14.*u1.*u3-M.*t3.*t7.*t9.*t10.*t67.*u1.*u2-I1.*M.*f2.*r.*t2.*t3.*t7.*t9.*t13-I3.*l.*r.*t7.*t8.*t9.*t13.*t14.*t68-I3.*l.*r.*t7.*t8.*t9.*t13.*t14.*t69-I2.*l.*r.*t3.*t6.*t8.*t10.*u2.*u3-I2.*l.*r.*t3.*t8.*t9.*t10.*u2.*u3-I2.*t2.*t3.*t6.*t8.*t9.*t13.*u2.*u3.*3.0-I3.*t3.*t6.*t7.*t8.*t9.*t10.*u1.*u2-l.*r.*t2.*t3.*t7.*t8.*t9.*t13.*tau3.*2.0-l.*r.*t2.*t8.*t9.*t10.*t13.*t14.*tau2.*2.0-l.*r.*t3.*t8.*t9.*t10.*t13.*t14.*tau1.*2.0-I2.*I3.*M.*l.*r.*t3.*t10.*u2.*u3-I2.*I3.*M.*t2.*t3.*t9.*t13.*u2.*u3-I1.*l.*r.*t2.*t8.*t9.*t10.*t13.*u1.*u3-I1.*l.*r.*t2.*t3.*t7.*t8.*t9.*t13.*u1.*u2-I3.*l.*r.*t2.*t3.*t7.*t8.*t9.*t13.*u1.*u2-I3.*l.*r.*t2.*t8.*t9.*t10.*t13.*t14.*u1.*u3.*2.0-I2.*l.*r.*t3.*t8.*t9.*t10.*t13.*t14.*u2.*u3.*2.0;
t169 = I1.*I3.*tau2p1;
t170 = I1.*M.*t6.*tau2p1;
t171 = I1.*t12.*u1p1.*u3p1;
t172 = I3.*M.*t9.*tau2p1;
t173 = t8.*t11.*t70.*tau2p1;
t174 = I1.*t8.*t29.*u1p1.*u3p1;
t175 = M.*t9.*t12.*u1p1.*u3p1;
t176 = I3.*M.*f1p1.*r.*t9.*t22;
t177 = t6.*t8.*t9.*t70.*tau2p1;
t178 = I1.*g.*l.*t6.*t8.*t22;
t179 = I1.*M.*t9.*t70.*tau2p1;
t180 = I1.*I3.*f1p1.*r.*t22;
t181 = I3.*M.*f3p1.*r.*t4.*t9.*t20.*t70;
t182 = I3.*t8.*t11.*t70.*u1p1.*u3p1;
t183 = I1.*l.*r.*t6.*t8.*t22.*t72;
t184 = I3.*l.*r.*t8.*t9.*t22.*t72;
t185 = I1.*l.*r.*t6.*t8.*t22.*t73;
t186 = I3.*l.*r.*t8.*t9.*t22.*t73;
t187 = I1.*M.*f1p1.*r.*t6.*t22;
t188 = l.*r.*t4.*t8.*t9.*t20.*t70.*tau2p1.*2.0;
t189 = I3.*M.*f3p1.*l.*t9.*t70;
t190 = I3.*t6.*t8.*t9.*u1p1.*u3p1;
t191 = I3.*g.*l.*t8.*t9.*t22;
t192 = I1.*I3.*M.*t6.*u1p1.*u3p1.*2.0;
t193 = I1.*I3.*M.*g.*l.*t22;
t194 = I1.*I3.*M.*t9.*t70.*u1p1.*u3p1;
t195 = I1.*t8.*t11.*t70.*t71.*u1p1.*u3p1;
t196 = I1.*M.*f1p1.*r.*t9.*t22.*t70.*t71;
t197 = I1.*I3.*M.*l.*r.*t22.*t72;
t198 = I1.*I3.*M.*l.*r.*t22.*t73;
t199 = I1.*M.*l.*r.*t4.*t20.*tau2p1.*2.0;
t200 = I1.*l.*r.*t4.*t6.*t8.*t20.*u1p1.*u3p1.*3.0;
t201 = I3.*l.*r.*t4.*t8.*t9.*t20.*u1p1.*u3p1;
t202 = I1.*M.*f1p1.*l.*t4.*t9.*t20.*t22.*2.0;
t203 = I3.*l.*r.*t5.*t8.*t9.*t20.*u2p1.*u3p1;
t204 = I3.*l.*r.*t5.*t6.*t8.*t20.*u2p1.*u3p1;
t205 = I3.*M.*f2p1.*l.*t5.*t9.*t20.*t22;
t206 = I1.*t6.*t8.*t9.*t70.*t71.*u1p1.*u3p1.*3.0;
t207 = I1.*g.*l.*t8.*t9.*t22.*t70.*t71;
t208 = I1.*t4.*t6.*t8.*t9.*t20.*t22.*t72.*2.0;
t209 = I1.*t4.*t6.*t8.*t9.*t20.*t22.*t73.*2.0;
t210 = I1.*I3.*M.*t9.*t70.*t71.*u1p1.*u3p1;
t211 = I3.*l.*r.*t4.*t8.*t9.*t20.*t70.*u1p1.*u3p1;
t212 = I1.*g.*r.*t4.*t6.*t8.*t20.*t22.*2.0;
t213 = I2.*t5.*t8.*t11.*t20.*t22.*u1p1.*u2p1;
t214 = I1.*l.*r.*t8.*t9.*t22.*t70.*t71.*t72;
t215 = I1.*l.*r.*t8.*t9.*t22.*t70.*t71.*t73;
t216 = M.*l.*r.*t5.*t12.*t20.*u2p1.*u3p1;
t217 = I1.*l.*r.*t4.*t8.*t9.*t20.*t70.*t71.*u1p1.*u3p1.*2.0;
t218 = M.*t4.*t5.*t9.*t12.*t70.*u2p1.*u3p1;
t219 = I3.*t4.*t5.*t8.*t11.*t70.*u2p1.*u3p1;
t220 = I3.*M.*f2p1.*r.*t4.*t5.*t9.*t22.*t70;
t221 = I2.*t5.*t6.*t8.*t9.*t20.*t22.*u1p1.*u2p1;
t222 = I1.*I2.*M.*t5.*t9.*t20.*t22.*u1p1.*u2p1;
t223 = I1.*I3.*M.*l.*r.*t4.*t20.*u1p1.*u3p1.*3.0;
t224 = I3.*t4.*t5.*t6.*t8.*t9.*t70.*u2p1.*u3p1.*3.0;
t225 = I3.*l.*r.*t5.*t8.*t9.*t20.*t70.*t71.*u2p1.*u3p1.*2.0;
t226 = I2.*l.*r.*t4.*t5.*t8.*t9.*t22.*t70.*u1p1.*u2p1.*2.0;
t227 = t169+t170+t171+t172+t173+t174+t175+t176+t177+t178+t179+t180+t181+t182+t183+t184+t185+t186+t187+t188+t189+t190+t191+t192+t193+t194+t195+t196+t197+t198+t199+t200+t201+t202+t203+t204+t205+t206+t207+t208+t209+t210+t211+t212+t213+t214+t215+t216+t217+t218+t219+t220+t221+t222+t223+t224+t225+t226-I3.*t67.*u1p1.*u3p1-I1.*I3.*f3p1.*l-I1.*M.*f3p1.*l.*t6-I3.*M.*f3p1.*l.*t9-M.*t6.*t67.*u1p1.*u3p1-t8.*t11.*t70.*t71.*tau2p1-t5.*t8.*t11.*t20.*t22.*tau3p1-t4.*t5.*t8.*t11.*t70.*tau1p1-t6.*t8.*t9.*t70.*t71.*tau2p1-I1.*I3.*M.*t9.*u1p1.*u3p1-I1.*M.*f3p1.*l.*t9.*t70-I1.*I3.*f3p1.*r.*t4.*t20-I3.*M.*t9.*t70.*t71.*tau2p1-I1.*t8.*t11.*t70.*u1p1.*u3p1-M.*t9.*t67.*t70.*u1p1.*u3p1-I1.*M.*f3p1.*l.*t9.*t70.*t71.*2.0-I1.*M.*f3p1.*r.*t4.*t6.*t20.*3.0-I3.*M.*f3p1.*r.*t4.*t9.*t20-I3.*M.*l.*r.*t5.*t20.*tau1p1-I1.*M.*t5.*t9.*t20.*t22.*tau3p1-I3.*M.*t4.*t5.*t9.*t70.*tau1p1-I3.*t8.*t11.*t70.*t71.*u1p1.*u3p1-M.*t9.*t12.*t70.*t71.*u1p1.*u3p1-l.*r.*t5.*t6.*t8.*t20.*tau1p1-l.*r.*t5.*t8.*t9.*t20.*tau1p1-t5.*t6.*t8.*t9.*t20.*t22.*tau3p1-t4.*t5.*t6.*t8.*t9.*t70.*tau1p1.*3.0-I1.*M.*f2p1.*l.*t5.*t9.*t20.*t22-I1.*M.*f3p1.*r.*t4.*t9.*t20.*t70-I3.*M.*f1p1.*r.*t9.*t22.*t70.*t71-I3.*g.*l.*t8.*t9.*t22.*t70.*t71-M.*l.*r.*t4.*t20.*t67.*u1p1.*u3p1.*2.0-I1.*t5.*t8.*t11.*t20.*t22.*u1p1.*u2p1-I2.*t4.*t5.*t8.*t11.*t70.*u2p1.*u3p1-I3.*t6.*t8.*t9.*t70.*t71.*u1p1.*u3p1-M.*t5.*t9.*t20.*t22.*t67.*u1p1.*u2p1-I1.*M.*f2p1.*r.*t4.*t5.*t9.*t22.*t70-I3.*l.*r.*t8.*t9.*t22.*t70.*t71.*t72-I3.*l.*r.*t8.*t9.*t22.*t70.*t71.*t73-I2.*l.*r.*t5.*t6.*t8.*t20.*u2p1.*u3p1-I2.*l.*r.*t5.*t8.*t9.*t20.*u2p1.*u3p1-I3.*t5.*t6.*t8.*t9.*t20.*t22.*u1p1.*u2p1-I2.*t4.*t5.*t6.*t8.*t9.*t70.*u2p1.*u3p1.*3.0-l.*r.*t4.*t5.*t8.*t9.*t22.*t70.*tau3p1.*2.0-l.*r.*t5.*t8.*t9.*t20.*t70.*t71.*tau1p1.*2.0-l.*r.*t4.*t8.*t9.*t20.*t70.*t71.*tau2p1.*2.0-I2.*I3.*M.*l.*r.*t5.*t20.*u2p1.*u3p1-I2.*I3.*M.*t4.*t5.*t9.*t70.*u2p1.*u3p1-I1.*l.*r.*t4.*t8.*t9.*t20.*t70.*u1p1.*u3p1-I1.*l.*r.*t4.*t5.*t8.*t9.*t22.*t70.*u1p1.*u2p1-I3.*l.*r.*t4.*t5.*t8.*t9.*t22.*t70.*u1p1.*u2p1-I3.*l.*r.*t4.*t8.*t9.*t20.*t70.*t71.*u1p1.*u3p1.*2.0-I2.*l.*r.*t5.*t8.*t9.*t20.*t70.*t71.*u2p1.*u3p1.*2.0;
t228 = t53.*t160-t89.*t227;
t229 = dt.*t228.*(1.0./8.0);
t230 = t100+t101+t229;
t231 = u1.*(1.0./2.0);
t232 = u1p1.*(1.0./2.0);
t233 = t8.*t29.*tau1;
t234 = I2.*I3.*tau1;
t235 = t6.*t8.*t9.*tau1;
t236 = I2.*M.*t6.*tau1;
t237 = I3.*M.*t6.*tau1;
t238 = I3.*t54.*u2.*u3;
t239 = I3.*M.*t9.*tau1;
t240 = M.*t6.*t54.*u2.*u3;
t241 = I2.*t8.*t29.*u2.*u3;
t242 = I2.*M.*t9.*t13.*tau1;
t243 = l.*r.*t7.*t8.*t9.*tau3;
t244 = l.*r.*t6.*t7.*t8.*tau3;
t245 = t8.*t11.*t13.*t14.*tau1;
t246 = l.*r.*t2.*t8.*t9.*t10.*t13.*t14.*tau1.*2.0;
t247 = I2.*M.*l.*r.*t7.*tau3;
t248 = M.*t9.*t13.*t54.*u2.*u3;
t249 = M.*t9.*t12.*t13.*u2.*u3;
t250 = I3.*M.*f2.*r.*t7.*t9.*t13;
t251 = I3.*M.*f3.*r.*t3.*t9.*t10.*t13;
t252 = l.*r.*t2.*t8.*t9.*t10.*tau1.*2.0;
t253 = l.*r.*t2.*t6.*t8.*t10.*tau1.*4.0;
t254 = t6.*t8.*t9.*t13.*t14.*tau1.*5.0;
t255 = I3.*M.*t9.*t13.*t14.*tau1;
t256 = I2.*t6.*t8.*t9.*u2.*u3;
t257 = t2.*t7.*t8.*t10.*t11.*tau3;
t258 = I2.*I3.*M.*t9.*u2.*u3;
t259 = l.*r.*t7.*t8.*t9.*t13.*t14.*tau3.*2.0;
t260 = I1.*l.*r.*t7.*t8.*t9.*u1.*u2;
t261 = I1.*l.*r.*t6.*t7.*t8.*u1.*u2;
t262 = I3.*l.*r.*t7.*t8.*t9.*u1.*u2;
t263 = I3.*l.*r.*t6.*t7.*t8.*u1.*u2;
t264 = I2.*t8.*t11.*t13.*t14.*u2.*u3;
t265 = I2.*M.*f2.*r.*t7.*t9.*t13.*t14;
t266 = t2.*t6.*t7.*t8.*t9.*t10.*tau3.*3.0;
t267 = I2.*M.*t2.*t7.*t9.*t10.*tau3;
t268 = I2.*M.*l.*r.*t2.*t10.*tau1.*2.0;
t269 = I3.*M.*l.*r.*t2.*t10.*tau1.*2.0;
t270 = I2.*l.*r.*t2.*t8.*t9.*t10.*u2.*u3.*2.0;
t271 = I2.*l.*r.*t2.*t6.*t8.*t10.*u2.*u3.*4.0;
t272 = I1.*l.*r.*t3.*t8.*t9.*t10.*u1.*u3;
t273 = I1.*l.*r.*t3.*t6.*t8.*t10.*u1.*u3;
t274 = I2.*l.*r.*t3.*t6.*t8.*t10.*u1.*u3;
t275 = I2.*M.*f1.*l.*t3.*t7.*t9.*t10;
t276 = I2.*t6.*t8.*t9.*t13.*t14.*u2.*u3.*5.0;
t277 = I2.*I3.*M.*t9.*t13.*t14.*u2.*u3;
t278 = I2.*t3.*t6.*t7.*t8.*t9.*t10.*t68;
t279 = I2.*t3.*t6.*t7.*t8.*t9.*t10.*t69;
t280 = I2.*l.*r.*t7.*t8.*t9.*t13.*u1.*u2;
t281 = I2.*l.*r.*t3.*t8.*t9.*t10.*t13.*u1.*u3;
t282 = I2.*g.*r.*t3.*t6.*t7.*t8.*t10;
t283 = I1.*t2.*t7.*t8.*t10.*t11.*u1.*u2;
t284 = M.*l.*r.*t2.*t10.*t54.*u2.*u3.*2.0;
t285 = I2.*l.*r.*t2.*t8.*t9.*t10.*t13.*t14.*u2.*u3.*2.0;
t286 = I1.*I2.*M.*l.*r.*t7.*u1.*u2;
t287 = I2.*I3.*M.*l.*r.*t7.*u1.*u2;
t288 = I1.*t2.*t3.*t8.*t11.*t13.*u1.*u3;
t289 = I2.*M.*f1.*r.*t2.*t3.*t7.*t9.*t13;
t290 = I2.*l.*r.*t2.*t3.*t7.*t8.*t9.*t13.*t68;
t291 = I2.*l.*r.*t2.*t3.*t7.*t8.*t9.*t13.*t69;
t292 = I1.*t2.*t6.*t7.*t8.*t9.*t10.*u1.*u2.*3.0;
t293 = I3.*t2.*t6.*t7.*t8.*t9.*t10.*u1.*u2.*2.0;
t294 = I1.*I2.*M.*t2.*t7.*t9.*t10.*u1.*u2;
t295 = I1.*I3.*M.*l.*r.*t3.*t10.*u1.*u3;
t296 = I2.*I3.*M.*l.*r.*t3.*t10.*u1.*u3;
t297 = I1.*l.*r.*t7.*t8.*t9.*t13.*t14.*u1.*u2.*2.0;
t298 = I1.*t2.*t3.*t6.*t8.*t9.*t13.*u1.*u3.*3.0;
t299 = I2.*t2.*t3.*t6.*t8.*t9.*t13.*u1.*u3.*2.0;
t300 = I3.*l.*r.*t7.*t8.*t9.*t13.*t14.*u1.*u2;
t301 = I1.*l.*r.*t3.*t8.*t9.*t10.*t13.*t14.*u1.*u3.*2.0;
t302 = I2.*g.*l.*t2.*t3.*t7.*t8.*t9.*t13;
t303 = I1.*I3.*M.*t2.*t3.*t9.*t13.*u1.*u3;
t304 = t233+t234+t235+t236+t237+t238+t239+t240+t241+t242+t243+t244+t245+t246+t247+t248+t249+t250+t251+t252+t253+t254+t255+t256+t257+t258+t259+t260+t261+t262+t263+t264+t265+t266+t267+t268+t269+t270+t271+t272+t273+t274+t275+t276+t277+t278+t279+t280+t281+t282+t283+t284+t285+t286+t287+t288+t289+t290+t291+t292+t293+t294+t295+t296+t297+t298+t299+t300+t301+t302+t303-I2.*t12.*u2.*u3-I2.*I3.*f2.*r.*t7-I3.*M.*t9.*t13.*tau1-I3.*t8.*t29.*u2.*u3-M.*t6.*t12.*u2.*u3-M.*t9.*t12.*u2.*u3-t2.*t3.*t8.*t11.*t13.*tau2-I2.*I3.*f3.*r.*t3.*t10-I3.*M.*f2.*r.*t6.*t7-I3.*M.*f2.*r.*t7.*t9-I3.*t6.*t8.*t9.*u2.*u3-I2.*I3.*M.*t9.*t13.*u2.*u3.*2.0-I2.*M.*f3.*r.*t3.*t6.*t10-I3.*M.*f3.*r.*t3.*t9.*t10-I2.*M.*f2.*r.*t7.*t9.*t13-I3.*M.*l.*r.*t3.*t10.*tau2-I3.*M.*t2.*t3.*t9.*t13.*tau2-M.*l.*r.*t7.*t54.*u1.*u2-I3.*t8.*t11.*t13.*t14.*u2.*u3-M.*t9.*t12.*t13.*t14.*u2.*u3-l.*r.*t3.*t6.*t8.*t10.*tau2-l.*r.*t3.*t8.*t9.*t10.*tau2-t2.*t3.*t6.*t8.*t9.*t13.*tau2.*3.0-I2.*M.*f3.*l.*t2.*t3.*t9.*t13.*2.0-I3.*M.*f1.*l.*t3.*t7.*t9.*t10-I3.*M.*f2.*l.*t2.*t7.*t9.*t10.*2.0-I2.*M.*f3.*r.*t3.*t9.*t10.*t13-I3.*M.*f2.*r.*t7.*t9.*t13.*t14-I3.*g.*r.*t3.*t6.*t7.*t8.*t10-I2.*l.*r.*t6.*t7.*t8.*u1.*u2-I2.*l.*r.*t7.*t8.*t9.*u1.*u2-M.*l.*r.*t2.*t10.*t12.*u2.*u3.*2.0-M.*l.*r.*t3.*t10.*t12.*u1.*u3-I3.*t3.*t6.*t7.*t8.*t9.*t10.*t68-I3.*t3.*t6.*t7.*t8.*t9.*t10.*t69-I2.*t2.*t7.*t8.*t10.*t11.*u1.*u2-I3.*t2.*t3.*t8.*t11.*t13.*u1.*u3-I3.*t6.*t8.*t9.*t13.*t14.*u2.*u3.*5.0-M.*t2.*t3.*t9.*t12.*t13.*u1.*u3-M.*t2.*t7.*t9.*t10.*t54.*u1.*u2-I3.*M.*f1.*r.*t2.*t3.*t7.*t9.*t13-I3.*g.*l.*t2.*t3.*t7.*t8.*t9.*t13-I3.*l.*r.*t2.*t6.*t8.*t10.*u2.*u3.*4.0-I3.*l.*r.*t3.*t6.*t8.*t10.*u1.*u3-I3.*l.*r.*t2.*t8.*t9.*t10.*u2.*u3.*2.0-I3.*l.*r.*t7.*t8.*t9.*t13.*u1.*u2-I2.*t2.*t6.*t7.*t8.*t9.*t10.*u1.*u2.*3.0-I3.*t2.*t3.*t6.*t8.*t9.*t13.*u1.*u3.*3.0-l.*r.*t3.*t8.*t9.*t10.*t13.*t14.*tau2.*2.0-I3.*l.*r.*t2.*t3.*t7.*t8.*t9.*t13.*t68-I3.*l.*r.*t2.*t3.*t7.*t8.*t9.*t13.*t69-I3.*l.*r.*t3.*t8.*t9.*t10.*t13.*u1.*u3-I2.*l.*r.*t7.*t8.*t9.*t13.*t14.*u1.*u2.*3.0-I3.*l.*r.*t2.*t8.*t9.*t10.*t13.*t14.*u2.*u3.*2.0-I3.*l.*r.*t3.*t8.*t9.*t10.*t13.*t14.*u1.*u3.*2.0;
t305 = t53.*t304;
t306 = t8.*t29.*tau1p1;
t307 = I2.*I3.*tau1p1;
t308 = t6.*t8.*t9.*tau1p1;
t309 = I2.*M.*t6.*tau1p1;
t310 = I3.*M.*t6.*tau1p1;
t311 = I3.*t54.*u2p1.*u3p1;
t312 = I3.*M.*t9.*tau1p1;
t313 = M.*t6.*t54.*u2p1.*u3p1;
t314 = I2.*t8.*t29.*u2p1.*u3p1;
t315 = I2.*M.*t9.*t70.*tau1p1;
t316 = l.*r.*t8.*t9.*t22.*tau3p1;
t317 = l.*r.*t6.*t8.*t22.*tau3p1;
t318 = t8.*t11.*t70.*t71.*tau1p1;
t319 = l.*r.*t4.*t8.*t9.*t20.*t70.*t71.*tau1p1.*2.0;
t320 = I2.*M.*l.*r.*t22.*tau3p1;
t321 = M.*t9.*t54.*t70.*u2p1.*u3p1;
t322 = M.*t9.*t12.*t70.*u2p1.*u3p1;
t323 = I3.*M.*f2p1.*r.*t9.*t22.*t70;
t324 = I3.*M.*f3p1.*r.*t5.*t9.*t20.*t70;
t325 = l.*r.*t4.*t8.*t9.*t20.*tau1p1.*2.0;
t326 = l.*r.*t4.*t6.*t8.*t20.*tau1p1.*4.0;
t327 = t6.*t8.*t9.*t70.*t71.*tau1p1.*5.0;
t328 = I3.*M.*t9.*t70.*t71.*tau1p1;
t329 = I2.*t6.*t8.*t9.*u2p1.*u3p1;
t330 = t4.*t8.*t11.*t20.*t22.*tau3p1;
t331 = I2.*I3.*M.*t9.*u2p1.*u3p1;
t332 = l.*r.*t8.*t9.*t22.*t70.*t71.*tau3p1.*2.0;
t333 = I1.*l.*r.*t8.*t9.*t22.*u1p1.*u2p1;
t334 = I1.*l.*r.*t6.*t8.*t22.*u1p1.*u2p1;
t335 = I3.*l.*r.*t8.*t9.*t22.*u1p1.*u2p1;
t336 = I3.*l.*r.*t6.*t8.*t22.*u1p1.*u2p1;
t337 = I2.*t8.*t11.*t70.*t71.*u2p1.*u3p1;
t338 = I2.*M.*f2p1.*r.*t9.*t22.*t70.*t71;
t339 = t4.*t6.*t8.*t9.*t20.*t22.*tau3p1.*3.0;
t340 = I2.*M.*t4.*t9.*t20.*t22.*tau3p1;
t341 = I2.*M.*l.*r.*t4.*t20.*tau1p1.*2.0;
t342 = I3.*M.*l.*r.*t4.*t20.*tau1p1.*2.0;
t343 = I2.*l.*r.*t4.*t8.*t9.*t20.*u2p1.*u3p1.*2.0;
t344 = I2.*l.*r.*t4.*t6.*t8.*t20.*u2p1.*u3p1.*4.0;
t345 = I1.*l.*r.*t5.*t8.*t9.*t20.*u1p1.*u3p1;
t346 = I1.*l.*r.*t5.*t6.*t8.*t20.*u1p1.*u3p1;
t347 = I2.*l.*r.*t5.*t6.*t8.*t20.*u1p1.*u3p1;
t348 = I2.*M.*f1p1.*l.*t5.*t9.*t20.*t22;
t349 = I2.*t6.*t8.*t9.*t70.*t71.*u2p1.*u3p1.*5.0;
t350 = I2.*I3.*M.*t9.*t70.*t71.*u2p1.*u3p1;
t351 = I2.*t5.*t6.*t8.*t9.*t20.*t22.*t72;
t352 = I2.*t5.*t6.*t8.*t9.*t20.*t22.*t73;
t353 = I2.*l.*r.*t8.*t9.*t22.*t70.*u1p1.*u2p1;
t354 = I2.*l.*r.*t5.*t8.*t9.*t20.*t70.*u1p1.*u3p1;
t355 = I2.*g.*r.*t5.*t6.*t8.*t20.*t22;
t356 = I1.*t4.*t8.*t11.*t20.*t22.*u1p1.*u2p1;
t357 = M.*l.*r.*t4.*t20.*t54.*u2p1.*u3p1.*2.0;
t358 = I2.*l.*r.*t4.*t8.*t9.*t20.*t70.*t71.*u2p1.*u3p1.*2.0;
t359 = I1.*I2.*M.*l.*r.*t22.*u1p1.*u2p1;
t360 = I2.*I3.*M.*l.*r.*t22.*u1p1.*u2p1;
t361 = I1.*t4.*t5.*t8.*t11.*t70.*u1p1.*u3p1;
t362 = I2.*M.*f1p1.*r.*t4.*t5.*t9.*t22.*t70;
t363 = I2.*l.*r.*t4.*t5.*t8.*t9.*t22.*t70.*t72;
t364 = I2.*l.*r.*t4.*t5.*t8.*t9.*t22.*t70.*t73;
t365 = I1.*t4.*t6.*t8.*t9.*t20.*t22.*u1p1.*u2p1.*3.0;
t366 = I3.*t4.*t6.*t8.*t9.*t20.*t22.*u1p1.*u2p1.*2.0;
t367 = I1.*I2.*M.*t4.*t9.*t20.*t22.*u1p1.*u2p1;
t368 = I1.*I3.*M.*l.*r.*t5.*t20.*u1p1.*u3p1;
t369 = I2.*I3.*M.*l.*r.*t5.*t20.*u1p1.*u3p1;
t370 = I1.*l.*r.*t8.*t9.*t22.*t70.*t71.*u1p1.*u2p1.*2.0;
t371 = I3.*l.*r.*t8.*t9.*t22.*t70.*t71.*u1p1.*u2p1;
t372 = I1.*t4.*t5.*t6.*t8.*t9.*t70.*u1p1.*u3p1.*3.0;
t373 = I2.*t4.*t5.*t6.*t8.*t9.*t70.*u1p1.*u3p1.*2.0;
t374 = I1.*l.*r.*t5.*t8.*t9.*t20.*t70.*t71.*u1p1.*u3p1.*2.0;
t375 = I2.*g.*l.*t4.*t5.*t8.*t9.*t22.*t70;
t376 = I1.*I3.*M.*t4.*t5.*t9.*t70.*u1p1.*u3p1;
t377 = t306+t307+t308+t309+t310+t311+t312+t313+t314+t315+t316+t317+t318+t319+t320+t321+t322+t323+t324+t325+t326+t327+t328+t329+t330+t331+t332+t333+t334+t335+t336+t337+t338+t339+t340+t341+t342+t343+t344+t345+t346+t347+t348+t349+t350+t351+t352+t353+t354+t355+t356+t357+t358+t359+t360+t361+t362+t363+t364+t365+t366+t367+t368+t369+t370+t371+t372+t373+t374+t375+t376-I2.*t12.*u2p1.*u3p1-I2.*I3.*f2p1.*r.*t22-I3.*M.*t9.*t70.*tau1p1-I3.*t8.*t29.*u2p1.*u3p1-M.*t6.*t12.*u2p1.*u3p1-M.*t9.*t12.*u2p1.*u3p1-t4.*t5.*t8.*t11.*t70.*tau2p1-I2.*I3.*f3p1.*r.*t5.*t20-I3.*M.*f2p1.*r.*t6.*t22-I3.*M.*f2p1.*r.*t9.*t22-I3.*t6.*t8.*t9.*u2p1.*u3p1-I2.*I3.*M.*t9.*t70.*u2p1.*u3p1.*2.0-I2.*M.*f3p1.*r.*t5.*t6.*t20-I3.*M.*f3p1.*r.*t5.*t9.*t20-I2.*M.*f2p1.*r.*t9.*t22.*t70-I3.*M.*l.*r.*t5.*t20.*tau2p1-I3.*M.*t4.*t5.*t9.*t70.*tau2p1-M.*l.*r.*t22.*t54.*u1p1.*u2p1-I3.*t8.*t11.*t70.*t71.*u2p1.*u3p1-M.*t9.*t12.*t70.*t71.*u2p1.*u3p1-l.*r.*t5.*t6.*t8.*t20.*tau2p1-l.*r.*t5.*t8.*t9.*t20.*tau2p1-t4.*t5.*t6.*t8.*t9.*t70.*tau2p1.*3.0-I3.*M.*f1p1.*l.*t5.*t9.*t20.*t22-I3.*M.*f2p1.*l.*t4.*t9.*t20.*t22.*2.0-I2.*M.*f3p1.*l.*t4.*t5.*t9.*t70.*2.0-I2.*M.*f3p1.*r.*t5.*t9.*t20.*t70-I3.*M.*f2p1.*r.*t9.*t22.*t70.*t71-I3.*g.*r.*t5.*t6.*t8.*t20.*t22-I2.*l.*r.*t6.*t8.*t22.*u1p1.*u2p1-I2.*l.*r.*t8.*t9.*t22.*u1p1.*u2p1-M.*l.*r.*t5.*t12.*t20.*u1p1.*u3p1-M.*l.*r.*t4.*t12.*t20.*u2p1.*u3p1.*2.0-I3.*t5.*t6.*t8.*t9.*t20.*t22.*t72-I3.*t5.*t6.*t8.*t9.*t20.*t22.*t73-I2.*t4.*t8.*t11.*t20.*t22.*u1p1.*u2p1-I3.*t4.*t5.*t8.*t11.*t70.*u1p1.*u3p1-I3.*t6.*t8.*t9.*t70.*t71.*u2p1.*u3p1.*5.0-M.*t4.*t9.*t20.*t22.*t54.*u1p1.*u2p1-M.*t4.*t5.*t9.*t12.*t70.*u1p1.*u3p1-I3.*M.*f1p1.*r.*t4.*t5.*t9.*t22.*t70-I3.*g.*l.*t4.*t5.*t8.*t9.*t22.*t70-I3.*l.*r.*t5.*t6.*t8.*t20.*u1p1.*u3p1-I3.*l.*r.*t4.*t6.*t8.*t20.*u2p1.*u3p1.*4.0-I3.*l.*r.*t4.*t8.*t9.*t20.*u2p1.*u3p1.*2.0-I3.*l.*r.*t8.*t9.*t22.*t70.*u1p1.*u2p1-I2.*t4.*t6.*t8.*t9.*t20.*t22.*u1p1.*u2p1.*3.0-I3.*t4.*t5.*t6.*t8.*t9.*t70.*u1p1.*u3p1.*3.0-l.*r.*t5.*t8.*t9.*t20.*t70.*t71.*tau2p1.*2.0-I3.*l.*r.*t4.*t5.*t8.*t9.*t22.*t70.*t72-I3.*l.*r.*t4.*t5.*t8.*t9.*t22.*t70.*t73-I3.*l.*r.*t5.*t8.*t9.*t20.*t70.*u1p1.*u3p1-I2.*l.*r.*t8.*t9.*t22.*t70.*t71.*u1p1.*u2p1.*3.0-I3.*l.*r.*t5.*t8.*t9.*t20.*t70.*t71.*u1p1.*u3p1.*2.0-I3.*l.*r.*t4.*t8.*t9.*t20.*t70.*t71.*u2p1.*u3p1.*2.0;
t378 = t305-t89.*t377;
t379 = dt.*t378.*(1.0./8.0);
t380 = t231+t232+t379;
out1 = (t74.*(dt.*t53.*(I1.*I3.*M.*l.*r.*t7.*u2.*2.0+I1.*l.*r.*t6.*t7.*t8.*u2.*2.0+I3.*l.*r.*t7.*t8.*t9.*u2.*2.0+M.*l.*r.*t3.*t10.*t12.*u3-I1.*t3.*t7.*t8.*t10.*t11.*u1-I2.*t2.*t3.*t8.*t11.*t13.*u3+I2.*t3.*t7.*t8.*t10.*t11.*u1+I3.*t2.*t3.*t8.*t11.*t13.*u3+M.*t2.*t3.*t9.*t12.*t13.*u3-M.*t3.*t7.*t9.*t10.*t67.*u1-I2.*I3.*M.*l.*r.*t3.*t10.*u3+I1.*I2.*M.*t3.*t7.*t9.*t10.*u1-I2.*I3.*M.*t2.*t3.*t9.*t13.*u3-I2.*l.*r.*t3.*t6.*t8.*t10.*u3+I3.*l.*r.*t3.*t6.*t8.*t10.*u3-I2.*l.*r.*t3.*t8.*t9.*t10.*u3+I3.*l.*r.*t3.*t8.*t9.*t10.*u3+I1.*t2.*t6.*t7.*t8.*t9.*t10.*u2.*4.0-I2.*t2.*t3.*t6.*t8.*t9.*t13.*u3.*3.0+I2.*t3.*t6.*t7.*t8.*t9.*t10.*u1+I3.*t2.*t3.*t6.*t8.*t9.*t13.*u3.*3.0-I3.*t3.*t6.*t7.*t8.*t9.*t10.*u1+I1.*l.*r.*t7.*t8.*t9.*t13.*t14.*u2.*2.0-I3.*l.*r.*t7.*t8.*t9.*t13.*t14.*u2.*2.0-I1.*l.*r.*t2.*t3.*t7.*t8.*t9.*t13.*u1+I2.*l.*r.*t2.*t3.*t7.*t8.*t9.*t13.*u1.*2.0-I3.*l.*r.*t2.*t3.*t7.*t8.*t9.*t13.*u1-I2.*l.*r.*t3.*t8.*t9.*t10.*t13.*t14.*u3.*2.0+I3.*l.*r.*t3.*t8.*t9.*t10.*t13.*t14.*u3.*2.0).*(1.0./8.0)+1.0./2.0)-dt.*t53.*t58.*(-I2.*t12.*u3+I3.*t54.*u3+I2.*t8.*t29.*u3-I3.*t8.*t29.*u3-M.*t6.*t12.*u3-M.*t9.*t12.*u3+M.*t6.*t54.*u3+I2.*I3.*M.*t9.*u3+I2.*t6.*t8.*t9.*u3-I3.*t6.*t8.*t9.*u3+M.*t9.*t12.*t13.*u3+M.*t9.*t13.*t54.*u3-I2.*I3.*M.*t9.*t13.*u3.*2.0-M.*l.*r.*t7.*t54.*u1+I2.*t8.*t11.*t13.*t14.*u3-I3.*t8.*t11.*t13.*t14.*u3-M.*t9.*t12.*t13.*t14.*u3+I1.*I2.*M.*l.*r.*t7.*u1+I2.*I3.*M.*l.*r.*t7.*u1+I2.*I3.*M.*t9.*t13.*t14.*u3+I1.*l.*r.*t6.*t7.*t8.*u1-I2.*l.*r.*t6.*t7.*t8.*u1+I3.*l.*r.*t6.*t7.*t8.*u1+I1.*l.*r.*t7.*t8.*t9.*u1-I2.*l.*r.*t7.*t8.*t9.*u1+I3.*l.*r.*t7.*t8.*t9.*u1-M.*l.*r.*t2.*t10.*t12.*u3.*2.0+M.*l.*r.*t2.*t10.*t54.*u3.*2.0+I1.*t2.*t7.*t8.*t10.*t11.*u1-I2.*t2.*t7.*t8.*t10.*t11.*u1+I2.*t6.*t8.*t9.*t13.*t14.*u3.*5.0-I3.*t6.*t8.*t9.*t13.*t14.*u3.*5.0-M.*t2.*t7.*t9.*t10.*t54.*u1+I1.*I2.*M.*t2.*t7.*t9.*t10.*u1+I2.*l.*r.*t2.*t6.*t8.*t10.*u3.*4.0-I3.*l.*r.*t2.*t6.*t8.*t10.*u3.*4.0+I2.*l.*r.*t2.*t8.*t9.*t10.*u3.*2.0-I3.*l.*r.*t2.*t8.*t9.*t10.*u3.*2.0+I2.*l.*r.*t7.*t8.*t9.*t13.*u1-I3.*l.*r.*t7.*t8.*t9.*t13.*u1+I1.*t2.*t6.*t7.*t8.*t9.*t10.*u1.*3.0-I2.*t2.*t6.*t7.*t8.*t9.*t10.*u1.*3.0+I3.*t2.*t6.*t7.*t8.*t9.*t10.*u1.*2.0+I2.*t3.*t6.*t7.*t8.*t9.*t10.*u2.*2.0-I3.*t3.*t6.*t7.*t8.*t9.*t10.*u2.*2.0+I1.*l.*r.*t7.*t8.*t9.*t13.*t14.*u1.*2.0-I2.*l.*r.*t7.*t8.*t9.*t13.*t14.*u1.*3.0+I3.*l.*r.*t7.*t8.*t9.*t13.*t14.*u1+I2.*l.*r.*t2.*t3.*t7.*t8.*t9.*t13.*u2.*2.0-I3.*l.*r.*t2.*t3.*t7.*t8.*t9.*t13.*u2.*2.0+I2.*l.*r.*t2.*t8.*t9.*t10.*t13.*t14.*u3.*2.0-I3.*l.*r.*t2.*t8.*t9.*t10.*t13.*t14.*u3.*2.0).*(1.0./8.0)+dt.*t3.*t7.*t17.*t58.*t230.*(1.0./8.0)+dt.*t3.*t7.*t17.*t74.*t380.*(1.0./8.0))./t97+t3.*t17.*(1.0./4.0)+dt.*t2.*1.0./t97.^2.*sin(t96).*(t74.*t230-t58.*t380).*(1.0./8.0);
