function out1 = dhi_dxp1_5_6(q1,q2,q3,q4,q5,u1,u2,u3,tau1,tau2,tau3,f1,f2,f3,q1p1,q2p1,q3p1,q4p1,q5p1,u1p1,u2p1,u3p1,tau1p1,tau2p1,tau3p1,f1p1,f2p1,f3p1,dt,M,l,r,g,I1,I2,I3)
%DHI_DXP1_5_6
%    OUT1 = DHI_DXP1_5_6(Q1,Q2,Q3,Q4,Q5,U1,U2,U3,TAU1,TAU2,TAU3,F1,F2,F3,Q1P1,Q2P1,Q3P1,Q4P1,Q5P1,U1P1,U2P1,U3P1,TAU1P1,TAU2P1,TAU3P1,F1P1,F2P1,F3P1,DT,M,L,R,G,I1,I2,I3)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    23-Nov-2018 21:05:46

t2 = l.^2;
t3 = M.^2;
t4 = r.^2;
t5 = t4.^2;
t6 = cos(q2);
t7 = t6.^2;
t8 = cos(q3);
t9 = t8.^2;
t10 = I1.^2;
t11 = t2.^2;
t12 = I3.^2;
t13 = sin(q2);
t14 = u2.^2;
t15 = u3.^2;
t16 = sin(q3);
t17 = I1.*t3.*t11;
t18 = I3.*t3.*t5;
t19 = I1.*I2.*I3;
t20 = I1.*I2.*M.*t2;
t21 = I1.*I3.*M.*t2;
t22 = I1.*I3.*M.*t4;
t23 = I2.*I3.*M.*t4;
t24 = cos(q2p1);
t25 = t24.^2;
t26 = I1.*t2.*t3.*t4;
t27 = I3.*t2.*t3.*t4;
t28 = cos(q3p1);
t29 = t28.^2;
t30 = sin(q2p1);
t31 = sin(q3p1);
t32 = u2p1.^2;
t33 = u3p1.^2;
t34 = 1.0./t6;
t35 = t8.*u1;
t42 = t16.*u2;
t36 = t35-t42;
t37 = 1.0./t24;
t38 = t28.*u1p1;
t44 = t31.*u2p1;
t39 = t38-t44;
t40 = q1.*(1.0./2.0);
t41 = q1p1.*(1.0./2.0);
t43 = t34.*t36;
t55 = t37.*t39;
t45 = t43-t55;
t46 = dt.*t45.*(1.0./8.0);
t47 = t40+t41+t46;
t48 = sin(t47);
t49 = q3.*(1.0./2.0);
t50 = q3p1.*(1.0./2.0);
t51 = t30.*t37.*t39;
t56 = t13.*t34.*t36;
t52 = t51-t56+u3-u3p1;
t53 = dt.*t52.*(1.0./8.0);
t54 = t49+t50+t53;
t57 = sin(t54);
t58 = cos(t47);
t59 = q2.*(1.0./2.0);
t60 = q2p1.*(1.0./2.0);
t61 = t8.*u2;
t62 = t16.*u1;
t67 = t28.*u2p1;
t68 = t31.*u1p1;
t63 = t61+t62-t67-t68;
t64 = dt.*t63.*(1.0./8.0);
t65 = t59+t60+t64;
t66 = cos(t54);
t69 = sin(t65);
t70 = I2.*t3.*t5.*t7;
t71 = I2.*t2.*t3.*t4.*t7;
t72 = I1.*I2.*M.*t4.*t7;
t73 = I1.*t3.*t5.*t7.*t9;
t74 = I1.*l.*r.*t3.*t4.*t6.*t8.*2.0;
t75 = I1.*l.*r.*t2.*t3.*t6.*t8.*4.0;
t76 = I3.*l.*r.*t3.*t4.*t6.*t8.*2.0;
t77 = I1.*t2.*t3.*t4.*t7.*t9.*5.0;
t78 = I1.*I3.*M.*t4.*t7.*t9;
t79 = I2.*l.*r.*t3.*t4.*t6.*t7.*t8.*2.0;
t80 = I1.*l.*r.*t3.*t4.*t6.*t7.*t8.*t9.*2.0;
t81 = I1.*I2.*M.*l.*r.*t6.*t8.*2.0;
t82 = I1.*I3.*M.*l.*r.*t6.*t8.*2.0;
t110 = I3.*t3.*t5.*t7;
t111 = I3.*t2.*t3.*t4.*t7;
t112 = I1.*I3.*M.*t4.*t7;
t113 = I2.*t3.*t5.*t7.*t9;
t114 = I2.*t2.*t3.*t4.*t7.*t9;
t115 = I2.*I3.*M.*t4.*t7.*t9;
t116 = I3.*l.*r.*t3.*t4.*t6.*t7.*t8.*2.0;
t117 = I2.*l.*r.*t3.*t4.*t6.*t7.*t8.*t9.*2.0;
t83 = t17+t18+t19+t20+t21+t22+t23+t26+t27+t70+t71+t72+t73+t74+t75+t76+t77+t78+t79+t80+t81+t82-t110-t111-t112-t113-t114-t115-t116-t117;
t84 = 1.0./t83;
t85 = I2.^2;
t86 = I2.*t3.*t5.*t25;
t87 = I2.*t2.*t3.*t4.*t25;
t88 = I1.*I2.*M.*t4.*t25;
t89 = I1.*t3.*t5.*t25.*t29;
t90 = I1.*l.*r.*t3.*t4.*t24.*t28.*2.0;
t91 = I1.*l.*r.*t2.*t3.*t24.*t28.*4.0;
t92 = I3.*l.*r.*t3.*t4.*t24.*t28.*2.0;
t93 = I1.*t2.*t3.*t4.*t25.*t29.*5.0;
t94 = I1.*I3.*M.*t4.*t25.*t29;
t95 = I2.*l.*r.*t3.*t4.*t24.*t25.*t28.*2.0;
t96 = I1.*l.*r.*t3.*t4.*t24.*t25.*t28.*t29.*2.0;
t97 = I1.*I2.*M.*l.*r.*t24.*t28.*2.0;
t98 = I1.*I3.*M.*l.*r.*t24.*t28.*2.0;
t102 = I3.*t3.*t5.*t25;
t103 = I3.*t2.*t3.*t4.*t25;
t104 = I1.*I3.*M.*t4.*t25;
t105 = I2.*t3.*t5.*t25.*t29;
t106 = I2.*t2.*t3.*t4.*t25.*t29;
t107 = I2.*I3.*M.*t4.*t25.*t29;
t108 = I3.*l.*r.*t3.*t4.*t24.*t25.*t28.*2.0;
t109 = I2.*l.*r.*t3.*t4.*t24.*t25.*t28.*t29.*2.0;
t99 = t17+t18+t19+t20+t21+t22+t23+t26+t27+t86+t87+t88+t89+t90+t91+t92+t93+t94+t95+t96+t97+t98-t102-t103-t104-t105-t106-t107-t108-t109;
t100 = 1.0./t99;
t101 = cos(t65);
t118 = u3.*(1.0./2.0);
t119 = u3p1.*(1.0./2.0);
t120 = t3.*t5.*tau3;
t121 = I1.*I2.*tau3;
t122 = t2.*t3.*t4.*tau3;
t123 = I1.*M.*f2.*l.*t2;
t124 = I1.*M.*t2.*tau3;
t125 = I2.*t10.*u1.*u2;
t126 = I1.*M.*t4.*tau3;
t127 = I2.*M.*t4.*tau3;
t128 = I1.*I2.*f2.*l;
t129 = M.*t2.*t10.*u1.*u2;
t130 = M.*t4.*t10.*u1.*u2;
t131 = I1.*t3.*t5.*u1.*u2;
t132 = I1.*M.*f2.*l.*t4;
t133 = l.*r.*t3.*t4.*t13.*tau1;
t134 = l.*r.*t2.*t3.*t13.*tau1;
t135 = I2.*M.*l.*r.*t13.*tau1;
t136 = I2.*M.*f2.*r.*t4.*t6.*t7.*t8;
t137 = I1.*g.*l.*t2.*t3.*t6.*t16;
t138 = I2.*t3.*t5.*t7.*u1.*u2;
t139 = l.*r.*t3.*t4.*t6.*t8.*tau3.*2.0;
t140 = I1.*I2.*f2.*r.*t6.*t8;
t141 = I1.*I2.*f1.*r.*t6.*t16;
t142 = I1.*M.*f2.*r.*t4.*t6.*t7.*t8.*t9;
t143 = I1.*M.*t4.*t7.*t9.*tau3;
t144 = I2.*M.*f2.*l.*t4.*t7;
t145 = I1.*M.*f2.*r.*t4.*t6.*t8;
t146 = t3.*t5.*t6.*t8.*t13.*tau1;
t147 = I2.*M.*f1.*r.*t4.*t6.*t16;
t148 = I1.*M.*f2.*l.*t4.*t7.*t9.*3.0;
t149 = M.*l.*r.*t13.*t85.*u2.*u3;
t150 = l.*r.*t3.*t4.*t7.*t9.*t13.*tau1.*2.0;
t151 = I1.*l.*r.*t2.*t3.*t6.*t14.*t16;
t152 = I1.*l.*r.*t2.*t3.*t6.*t15.*t16;
t153 = I2.*l.*r.*t3.*t4.*t6.*t14.*t16;
t154 = I2.*l.*r.*t3.*t4.*t6.*t15.*t16;
t155 = I1.*I2.*M.*t4.*t7.*u1.*u2;
t156 = I1.*M.*f2.*r.*t2.*t6.*t8.*3.0;
t157 = I1.*M.*f1.*r.*t2.*t6.*t16;
t158 = I2.*l.*r.*t3.*t4.*t13.*u2.*u3;
t159 = I2.*l.*r.*t2.*t3.*t13.*u2.*u3;
t160 = M.*t4.*t7.*t9.*t10.*u1.*u2;
t161 = M.*t4.*t7.*t9.*t85.*u1.*u2;
t162 = I1.*M.*f1.*r.*t4.*t6.*t7.*t9.*t16;
t163 = t2.*t3.*t4.*t6.*t8.*t13.*tau1.*3.0;
t164 = I2.*g.*l.*t3.*t4.*t6.*t16;
t165 = I2.*M.*t4.*t6.*t8.*t13.*tau1;
t166 = I1.*I2.*M.*g.*l.*t6.*t16;
t167 = I1.*M.*l.*r.*t6.*t8.*tau3.*2.0;
t168 = I1.*l.*r.*t3.*t4.*t6.*t8.*u1.*u2;
t169 = I1.*M.*f3.*l.*t4.*t6.*t13.*t16;
t170 = I2.*t2.*t3.*t4.*t7.*t9.*u1.*u2;
t171 = I1.*g.*l.*t3.*t4.*t6.*t7.*t9.*t16;
t172 = I2.*l.*r.*t3.*t4.*t6.*t7.*t8.*u1.*u2;
t173 = I1.*M.*f1.*l.*t4.*t7.*t8.*t16.*2.0;
t174 = M.*t4.*t6.*t8.*t13.*t85.*u2.*u3;
t175 = I2.*t3.*t5.*t6.*t8.*t13.*u2.*u3;
t176 = I1.*I2.*M.*l.*r.*t6.*t14.*t16;
t177 = I1.*I2.*M.*l.*r.*t6.*t15.*t16;
t178 = M.*t4.*t6.*t10.*t13.*t16.*u1.*u3;
t179 = I1.*t3.*t5.*t6.*t13.*t16.*u1.*u3;
t180 = I1.*t2.*t3.*t4.*t7.*t8.*t14.*t16.*2.0;
t181 = I1.*t2.*t3.*t4.*t7.*t8.*t15.*t16.*2.0;
t182 = I1.*l.*r.*t3.*t4.*t6.*t7.*t9.*t14.*t16;
t183 = I1.*l.*r.*t3.*t4.*t6.*t7.*t9.*t15.*t16;
t184 = M.*l.*r.*t6.*t8.*t10.*u1.*u2.*2.0;
t185 = I2.*l.*r.*t3.*t4.*t6.*t7.*t8.*t9.*u1.*u2;
t186 = I1.*g.*r.*t2.*t3.*t7.*t8.*t16.*2.0;
t187 = I1.*M.*f3.*r.*t4.*t7.*t8.*t13.*t16;
t188 = I2.*t2.*t3.*t4.*t6.*t8.*t13.*u2.*u3.*3.0;
t189 = I2.*t2.*t3.*t4.*t6.*t13.*t16.*u1.*u3;
t190 = I2.*l.*r.*t3.*t4.*t7.*t9.*t13.*u2.*u3.*2.0;
t191 = I1.*l.*r.*t3.*t4.*t7.*t8.*t13.*t16.*u1.*u3;
t192 = I2.*l.*r.*t3.*t4.*t7.*t8.*t13.*t16.*u1.*u3;
t193 = t120+t121+t122+t123+t124+t125+t126+t127+t128+t129+t130+t131+t132+t133+t134+t135+t136+t137+t138+t139+t140+t141+t142+t143+t144+t145+t146+t147+t148+t149+t150+t151+t152+t153+t154+t155+t156+t157+t158+t159+t160+t161+t162+t163+t164+t165+t166+t167+t168+t169+t170+t171+t172+t173+t174+t175+t176+t177+t178+t179+t180+t181+t182+t183+t184+t185+t186+t187+t188+t189+t190+t191+t192-I1.*t85.*u1.*u2-t3.*t5.*t7.*tau3-I1.*M.*t4.*t7.*tau3-I2.*t3.*t5.*u1.*u2-I1.*t3.*t11.*u1.*u2-M.*t4.*t85.*u1.*u2-t2.*t3.*t4.*t7.*tau3-t3.*t5.*t6.*t13.*t16.*tau2-I1.*I2.*M.*t2.*u1.*u2.*2.0-I1.*M.*f2.*l.*t4.*t7-I2.*M.*t4.*t7.*t9.*tau3-I2.*t2.*t3.*t4.*u1.*u2-I1.*t3.*t5.*t7.*u1.*u2-M.*t4.*t7.*t10.*u1.*u2-I2.*M.*f2.*l.*t4.*t7.*t9-I1.*M.*t4.*t6.*t13.*t16.*tau2-t2.*t3.*t4.*t6.*t13.*t16.*tau2-I2.*I3.*M.*l.*r.*t13.*u2.*u3-I1.*I2.*M.*t4.*t7.*t9.*u1.*u2.*2.0-I2.*M.*f3.*l.*t4.*t6.*t13.*t16-I1.*M.*f2.*r.*t4.*t6.*t7.*t8-I3.*l.*r.*t2.*t3.*t13.*u2.*u3-I3.*l.*r.*t3.*t4.*t13.*u2.*u3-I1.*t2.*t3.*t4.*t7.*t9.*u1.*u2.*3.0-I3.*t3.*t5.*t6.*t8.*t13.*u2.*u3-I3.*t3.*t5.*t6.*t13.*t16.*u1.*u3-l.*r.*t3.*t4.*t6.*t7.*t8.*tau3.*2.0-I2.*M.*f2.*r.*t4.*t6.*t7.*t8.*t9-I2.*M.*f1.*r.*t4.*t6.*t7.*t9.*t16-I2.*M.*f3.*r.*t4.*t7.*t8.*t13.*t16-I2.*g.*l.*t3.*t4.*t6.*t7.*t9.*t16-I1.*l.*r.*t2.*t3.*t6.*t8.*u1.*u2.*3.0-I2.*l.*r.*t3.*t4.*t6.*t8.*u1.*u2.*2.0-I3.*t2.*t3.*t4.*t6.*t8.*t13.*u2.*u3.*3.0-I3.*t2.*t3.*t4.*t6.*t13.*t16.*u1.*u3-l.*r.*t3.*t4.*t7.*t8.*t13.*t16.*tau2.*2.0-I1.*I2.*M.*l.*r.*t6.*t8.*u1.*u2.*3.0-I2.*I3.*M.*t4.*t6.*t8.*t13.*u2.*u3-I1.*I3.*M.*t4.*t6.*t13.*t16.*u1.*u3-I2.*l.*r.*t3.*t4.*t6.*t7.*t9.*t14.*t16-I2.*l.*r.*t3.*t4.*t6.*t7.*t9.*t15.*t16-I1.*l.*r.*t3.*t4.*t6.*t7.*t8.*u1.*u2-I3.*l.*r.*t3.*t4.*t7.*t9.*t13.*u2.*u3.*2.0-I1.*l.*r.*t3.*t4.*t6.*t7.*t8.*t9.*u1.*u2-I3.*l.*r.*t3.*t4.*t7.*t8.*t13.*t16.*u1.*u3.*2.0;
t194 = t84.*t193;
t195 = t3.*t5.*tau3p1;
t196 = I1.*I2.*tau3p1;
t197 = t2.*t3.*t4.*tau3p1;
t198 = I1.*M.*f2p1.*l.*t2;
t199 = I1.*M.*t2.*tau3p1;
t200 = I2.*t10.*u1p1.*u2p1;
t201 = I1.*M.*t4.*tau3p1;
t202 = I2.*M.*t4.*tau3p1;
t203 = I1.*I2.*f2p1.*l;
t204 = M.*t2.*t10.*u1p1.*u2p1;
t205 = M.*t4.*t10.*u1p1.*u2p1;
t206 = I1.*t3.*t5.*u1p1.*u2p1;
t207 = I1.*M.*f2p1.*l.*t4;
t208 = l.*r.*t3.*t4.*t30.*tau1p1;
t209 = l.*r.*t2.*t3.*t30.*tau1p1;
t210 = I2.*M.*l.*r.*t30.*tau1p1;
t211 = I2.*M.*f2p1.*r.*t4.*t24.*t25.*t28;
t212 = I1.*g.*l.*t2.*t3.*t24.*t31;
t213 = I2.*t3.*t5.*t25.*u1p1.*u2p1;
t214 = l.*r.*t3.*t4.*t24.*t28.*tau3p1.*2.0;
t215 = I1.*I2.*f2p1.*r.*t24.*t28;
t216 = I1.*I2.*f1p1.*r.*t24.*t31;
t217 = I1.*M.*f2p1.*r.*t4.*t24.*t25.*t28.*t29;
t218 = I1.*M.*t4.*t25.*t29.*tau3p1;
t219 = I2.*M.*f2p1.*l.*t4.*t25;
t220 = I1.*M.*f2p1.*r.*t4.*t24.*t28;
t221 = t3.*t5.*t24.*t28.*t30.*tau1p1;
t222 = I2.*M.*f1p1.*r.*t4.*t24.*t31;
t223 = I1.*M.*f2p1.*l.*t4.*t25.*t29.*3.0;
t224 = M.*l.*r.*t30.*t85.*u2p1.*u3p1;
t225 = l.*r.*t3.*t4.*t25.*t29.*t30.*tau1p1.*2.0;
t226 = I1.*l.*r.*t2.*t3.*t24.*t31.*t32;
t227 = I2.*l.*r.*t3.*t4.*t24.*t31.*t32;
t228 = I1.*l.*r.*t2.*t3.*t24.*t31.*t33;
t229 = I2.*l.*r.*t3.*t4.*t24.*t31.*t33;
t230 = I1.*I2.*M.*t4.*t25.*u1p1.*u2p1;
t231 = I1.*M.*f2p1.*r.*t2.*t24.*t28.*3.0;
t232 = I1.*M.*f1p1.*r.*t2.*t24.*t31;
t233 = I2.*l.*r.*t3.*t4.*t30.*u2p1.*u3p1;
t234 = I2.*l.*r.*t2.*t3.*t30.*u2p1.*u3p1;
t235 = M.*t4.*t10.*t25.*t29.*u1p1.*u2p1;
t236 = M.*t4.*t25.*t29.*t85.*u1p1.*u2p1;
t237 = I1.*M.*f1p1.*r.*t4.*t24.*t25.*t29.*t31;
t238 = t2.*t3.*t4.*t24.*t28.*t30.*tau1p1.*3.0;
t239 = I2.*g.*l.*t3.*t4.*t24.*t31;
t240 = I2.*M.*t4.*t24.*t28.*t30.*tau1p1;
t241 = I1.*I2.*M.*g.*l.*t24.*t31;
t242 = I1.*M.*l.*r.*t24.*t28.*tau3p1.*2.0;
t243 = I1.*l.*r.*t3.*t4.*t24.*t28.*u1p1.*u2p1;
t244 = I1.*M.*f3p1.*l.*t4.*t24.*t30.*t31;
t245 = I2.*t2.*t3.*t4.*t25.*t29.*u1p1.*u2p1;
t246 = I1.*g.*l.*t3.*t4.*t24.*t25.*t29.*t31;
t247 = I2.*l.*r.*t3.*t4.*t24.*t25.*t28.*u1p1.*u2p1;
t248 = I1.*M.*f1p1.*l.*t4.*t25.*t28.*t31.*2.0;
t249 = M.*t4.*t24.*t28.*t30.*t85.*u2p1.*u3p1;
t250 = I2.*t3.*t5.*t24.*t28.*t30.*u2p1.*u3p1;
t251 = I1.*I2.*M.*l.*r.*t24.*t31.*t32;
t252 = I1.*I2.*M.*l.*r.*t24.*t31.*t33;
t253 = M.*t4.*t10.*t24.*t30.*t31.*u1p1.*u3p1;
t254 = I1.*t3.*t5.*t24.*t30.*t31.*u1p1.*u3p1;
t255 = I1.*t2.*t3.*t4.*t25.*t28.*t31.*t32.*2.0;
t256 = I1.*l.*r.*t3.*t4.*t24.*t25.*t29.*t31.*t32;
t257 = I1.*t2.*t3.*t4.*t25.*t28.*t31.*t33.*2.0;
t258 = I1.*l.*r.*t3.*t4.*t24.*t25.*t29.*t31.*t33;
t259 = M.*l.*r.*t10.*t24.*t28.*u1p1.*u2p1.*2.0;
t260 = I2.*l.*r.*t3.*t4.*t24.*t25.*t28.*t29.*u1p1.*u2p1;
t261 = I1.*g.*r.*t2.*t3.*t25.*t28.*t31.*2.0;
t262 = I1.*M.*f3p1.*r.*t4.*t25.*t28.*t30.*t31;
t263 = I2.*t2.*t3.*t4.*t24.*t28.*t30.*u2p1.*u3p1.*3.0;
t264 = I2.*t2.*t3.*t4.*t24.*t30.*t31.*u1p1.*u3p1;
t265 = I2.*l.*r.*t3.*t4.*t25.*t29.*t30.*u2p1.*u3p1.*2.0;
t266 = I1.*l.*r.*t3.*t4.*t25.*t28.*t30.*t31.*u1p1.*u3p1;
t267 = I2.*l.*r.*t3.*t4.*t25.*t28.*t30.*t31.*u1p1.*u3p1;
t268 = t195+t196+t197+t198+t199+t200+t201+t202+t203+t204+t205+t206+t207+t208+t209+t210+t211+t212+t213+t214+t215+t216+t217+t218+t219+t220+t221+t222+t223+t224+t225+t226+t227+t228+t229+t230+t231+t232+t233+t234+t235+t236+t237+t238+t239+t240+t241+t242+t243+t244+t245+t246+t247+t248+t249+t250+t251+t252+t253+t254+t255+t256+t257+t258+t259+t260+t261+t262+t263+t264+t265+t266+t267-I1.*t85.*u1p1.*u2p1-t3.*t5.*t25.*tau3p1-I1.*M.*t4.*t25.*tau3p1-I2.*t3.*t5.*u1p1.*u2p1-I1.*t3.*t11.*u1p1.*u2p1-M.*t4.*t85.*u1p1.*u2p1-t2.*t3.*t4.*t25.*tau3p1-t3.*t5.*t24.*t30.*t31.*tau2p1-I1.*I2.*M.*t2.*u1p1.*u2p1.*2.0-I1.*M.*f2p1.*l.*t4.*t25-I2.*M.*t4.*t25.*t29.*tau3p1-I2.*t2.*t3.*t4.*u1p1.*u2p1-I1.*t3.*t5.*t25.*u1p1.*u2p1-M.*t4.*t10.*t25.*u1p1.*u2p1-I2.*M.*f2p1.*l.*t4.*t25.*t29-I1.*M.*t4.*t24.*t30.*t31.*tau2p1-t2.*t3.*t4.*t24.*t30.*t31.*tau2p1-I2.*I3.*M.*l.*r.*t30.*u2p1.*u3p1-I1.*I2.*M.*t4.*t25.*t29.*u1p1.*u2p1.*2.0-I2.*M.*f3p1.*l.*t4.*t24.*t30.*t31-I1.*M.*f2p1.*r.*t4.*t24.*t25.*t28-I3.*l.*r.*t2.*t3.*t30.*u2p1.*u3p1-I3.*l.*r.*t3.*t4.*t30.*u2p1.*u3p1-I1.*t2.*t3.*t4.*t25.*t29.*u1p1.*u2p1.*3.0-I3.*t3.*t5.*t24.*t30.*t31.*u1p1.*u3p1-I3.*t3.*t5.*t24.*t28.*t30.*u2p1.*u3p1-l.*r.*t3.*t4.*t24.*t25.*t28.*tau3p1.*2.0-I2.*M.*f1p1.*r.*t4.*t24.*t25.*t29.*t31-I2.*M.*f2p1.*r.*t4.*t24.*t25.*t28.*t29-I2.*M.*f3p1.*r.*t4.*t25.*t28.*t30.*t31-I2.*g.*l.*t3.*t4.*t24.*t25.*t29.*t31-I1.*l.*r.*t2.*t3.*t24.*t28.*u1p1.*u2p1.*3.0-I2.*l.*r.*t3.*t4.*t24.*t28.*u1p1.*u2p1.*2.0-I3.*t2.*t3.*t4.*t24.*t30.*t31.*u1p1.*u3p1-I3.*t2.*t3.*t4.*t24.*t28.*t30.*u2p1.*u3p1.*3.0-l.*r.*t3.*t4.*t25.*t28.*t30.*t31.*tau2p1.*2.0-I1.*I2.*M.*l.*r.*t24.*t28.*u1p1.*u2p1.*3.0-I1.*I3.*M.*t4.*t24.*t30.*t31.*u1p1.*u3p1-I2.*I3.*M.*t4.*t24.*t28.*t30.*u2p1.*u3p1-I2.*l.*r.*t3.*t4.*t24.*t25.*t29.*t31.*t32-I2.*l.*r.*t3.*t4.*t24.*t25.*t29.*t31.*t33-I1.*l.*r.*t3.*t4.*t24.*t25.*t28.*u1p1.*u2p1-I3.*l.*r.*t3.*t4.*t25.*t29.*t30.*u2p1.*u3p1.*2.0-I1.*l.*r.*t3.*t4.*t24.*t25.*t28.*t29.*u1p1.*u2p1-I3.*l.*r.*t3.*t4.*t25.*t28.*t30.*t31.*u1p1.*u3p1.*2.0;
t269 = t194-t100.*t268;
t270 = dt.*t269.*(1.0./8.0);
t271 = t118+t119+t270;
out1 = r.*(t31.*cos(q1p1)+t28.*t30.*sin(q1p1)).*(1.0./4.0)-r.*((u1.*(1.0./2.0)+u1p1.*(1.0./2.0)+dt.*(t84.*(I2.*I3.*tau1+t3.*t11.*tau1-I2.*t12.*u2.*u3+I3.*t85.*u2.*u3+t2.*t3.*t4.*tau1+I2.*M.*t2.*tau1+I3.*M.*t2.*tau1+I3.*M.*t4.*tau1-I2.*I3.*f2.*r.*t13+I2.*M.*t4.*t7.*tau1-I3.*M.*t4.*t7.*tau1+I2.*t3.*t11.*u2.*u3-I3.*t3.*t11.*u2.*u3-M.*t2.*t12.*u2.*u3-M.*t4.*t12.*u2.*u3+M.*t2.*t85.*u2.*u3+t3.*t5.*t7.*t9.*tau1+l.*r.*t2.*t3.*t13.*tau3+l.*r.*t3.*t4.*t13.*tau3+t2.*t3.*t4.*t7.*t9.*tau1.*5.0+t3.*t5.*t6.*t8.*t13.*tau3-t3.*t5.*t7.*t8.*t16.*tau2+I2.*I3.*M.*t4.*u2.*u3-I2.*I3.*f3.*r.*t6.*t16-I3.*M.*f2.*r.*t2.*t13-I3.*M.*f2.*r.*t4.*t13+I2.*M.*l.*r.*t13.*tau3+I3.*M.*t4.*t7.*t9.*tau1+I2.*t2.*t3.*t4.*u2.*u3-I3.*t2.*t3.*t4.*u2.*u3+M.*t4.*t7.*t12.*u2.*u3+M.*t4.*t7.*t85.*u2.*u3-I2.*I3.*M.*t4.*t7.*u2.*u3.*2.0-I2.*M.*f2.*r.*t4.*t7.*t13-I2.*M.*f3.*r.*t2.*t6.*t16+I3.*M.*f2.*r.*t4.*t7.*t13-I3.*M.*f3.*r.*t4.*t6.*t16+I2.*M.*l.*r.*t6.*t8.*tau1.*2.0+I3.*M.*l.*r.*t6.*t8.*tau1.*2.0-I3.*M.*l.*r.*t6.*t16.*tau2+I2.*M.*t4.*t6.*t8.*t13.*tau3-I3.*M.*t4.*t7.*t8.*t16.*tau2-M.*l.*r.*t13.*t85.*u1.*u2+I2.*t3.*t5.*t7.*t9.*u2.*u3-I3.*t3.*t5.*t7.*t9.*u2.*u3-M.*t4.*t7.*t9.*t12.*u2.*u3+l.*r.*t2.*t3.*t6.*t8.*tau1.*4.0+l.*r.*t3.*t4.*t6.*t8.*tau1.*2.0-l.*r.*t2.*t3.*t6.*t16.*tau2-l.*r.*t3.*t4.*t6.*t16.*tau2+t2.*t3.*t4.*t6.*t8.*t13.*tau3.*3.0-t2.*t3.*t4.*t7.*t8.*t16.*tau2.*3.0+I1.*I2.*M.*l.*r.*t13.*u1.*u2+I2.*I3.*M.*l.*r.*t13.*u1.*u2+I2.*I3.*M.*t4.*t7.*t9.*u2.*u3-I3.*M.*f2.*l.*t4.*t6.*t8.*t13.*2.0-I2.*M.*f3.*l.*t4.*t7.*t8.*t16.*2.0+I2.*M.*f1.*l.*t4.*t6.*t13.*t16-I3.*M.*f1.*l.*t4.*t6.*t13.*t16+I2.*M.*f2.*r.*t4.*t7.*t9.*t13-I2.*M.*f3.*r.*t4.*t6.*t7.*t16-I3.*M.*f2.*r.*t4.*t7.*t9.*t13+I3.*M.*f3.*r.*t4.*t6.*t7.*t16+I2.*g.*r.*t2.*t3.*t6.*t13.*t16-I3.*g.*r.*t2.*t3.*t6.*t13.*t16+I1.*l.*r.*t2.*t3.*t13.*u1.*u2-I2.*l.*r.*t2.*t3.*t13.*u1.*u2+I1.*l.*r.*t3.*t4.*t13.*u1.*u2+I3.*l.*r.*t2.*t3.*t13.*u1.*u2-I2.*l.*r.*t3.*t4.*t13.*u1.*u2+I3.*l.*r.*t3.*t4.*t13.*u1.*u2-M.*l.*r.*t6.*t8.*t12.*u2.*u3.*2.0-M.*l.*r.*t6.*t12.*t16.*u1.*u3+M.*l.*r.*t6.*t8.*t85.*u2.*u3.*2.0+I2.*t2.*t3.*t4.*t6.*t13.*t14.*t16+I2.*t2.*t3.*t4.*t6.*t13.*t15.*t16-I3.*t2.*t3.*t4.*t6.*t13.*t14.*t16-I3.*t2.*t3.*t4.*t6.*t13.*t15.*t16+I2.*t2.*t3.*t4.*t7.*t9.*u2.*u3.*5.0-I3.*t2.*t3.*t4.*t7.*t9.*u2.*u3.*5.0+I1.*t3.*t5.*t6.*t8.*t13.*u1.*u2-I2.*t3.*t5.*t6.*t8.*t13.*u1.*u2+I1.*t3.*t5.*t7.*t8.*t16.*u1.*u3-I3.*t3.*t5.*t7.*t8.*t16.*u1.*u3-M.*t4.*t7.*t8.*t12.*t16.*u1.*u3-M.*t4.*t6.*t8.*t13.*t85.*u1.*u2+l.*r.*t3.*t4.*t7.*t9.*t13.*tau3.*2.0+I2.*M.*f1.*r.*t4.*t7.*t8.*t13.*t16-I3.*M.*f1.*r.*t4.*t7.*t8.*t13.*t16+I2.*g.*l.*t3.*t4.*t7.*t8.*t13.*t16-I3.*g.*l.*t3.*t4.*t7.*t8.*t13.*t16+I2.*l.*r.*t2.*t3.*t6.*t8.*u2.*u3.*4.0-I3.*l.*r.*t2.*t3.*t6.*t8.*u2.*u3.*4.0+I2.*l.*r.*t3.*t4.*t6.*t8.*u2.*u3.*2.0-I3.*l.*r.*t3.*t4.*t6.*t8.*u2.*u3.*2.0+I1.*l.*r.*t2.*t3.*t6.*t16.*u1.*u3+I2.*l.*r.*t3.*t4.*t7.*t13.*u1.*u2+I2.*l.*r.*t2.*t3.*t6.*t16.*u1.*u3-I3.*l.*r.*t3.*t4.*t7.*t13.*u1.*u2+I1.*l.*r.*t3.*t4.*t6.*t16.*u1.*u3-I3.*l.*r.*t2.*t3.*t6.*t16.*u1.*u3+I1.*t2.*t3.*t4.*t6.*t8.*t13.*u1.*u2.*3.0-I2.*t2.*t3.*t4.*t6.*t8.*t13.*u1.*u2.*3.0+I3.*t2.*t3.*t4.*t6.*t8.*t13.*u1.*u2.*2.0+I1.*t2.*t3.*t4.*t7.*t8.*t16.*u1.*u3.*3.0+I2.*t2.*t3.*t4.*t7.*t8.*t16.*u1.*u3.*2.0-I3.*t2.*t3.*t4.*t7.*t8.*t16.*u1.*u3.*3.0+l.*r.*t3.*t4.*t6.*t7.*t8.*t9.*tau1.*2.0-l.*r.*t3.*t4.*t6.*t7.*t9.*t16.*tau2.*2.0+I1.*I3.*M.*l.*r.*t6.*t16.*u1.*u3+I2.*I3.*M.*l.*r.*t6.*t16.*u1.*u3+I1.*I2.*M.*t4.*t6.*t8.*t13.*u1.*u2+I1.*I3.*M.*t4.*t7.*t8.*t16.*u1.*u3+I2.*l.*r.*t3.*t4.*t7.*t8.*t13.*t14.*t16+I2.*l.*r.*t3.*t4.*t7.*t8.*t13.*t15.*t16-I3.*l.*r.*t3.*t4.*t7.*t8.*t13.*t14.*t16-I3.*l.*r.*t3.*t4.*t7.*t8.*t13.*t15.*t16+I1.*l.*r.*t3.*t4.*t7.*t9.*t13.*u1.*u2.*2.0-I2.*l.*r.*t3.*t4.*t7.*t9.*t13.*u1.*u2.*3.0+I2.*l.*r.*t3.*t4.*t6.*t7.*t16.*u1.*u3+I3.*l.*r.*t3.*t4.*t7.*t9.*t13.*u1.*u2-I3.*l.*r.*t3.*t4.*t6.*t7.*t16.*u1.*u3+I2.*l.*r.*t3.*t4.*t6.*t7.*t8.*t9.*u2.*u3.*2.0-I3.*l.*r.*t3.*t4.*t6.*t7.*t8.*t9.*u2.*u3.*2.0+I1.*l.*r.*t3.*t4.*t6.*t7.*t9.*t16.*u1.*u3.*2.0-I3.*l.*r.*t3.*t4.*t6.*t7.*t9.*t16.*u1.*u3.*2.0)-t100.*(I2.*I3.*tau1p1+t3.*t11.*tau1p1-I2.*t12.*u2p1.*u3p1+I3.*t85.*u2p1.*u3p1+t2.*t3.*t4.*tau1p1+I2.*M.*t2.*tau1p1+I3.*M.*t2.*tau1p1+I3.*M.*t4.*tau1p1-I2.*I3.*f2p1.*r.*t30+I2.*M.*t4.*t25.*tau1p1-I3.*M.*t4.*t25.*tau1p1+I2.*t3.*t11.*u2p1.*u3p1-I3.*t3.*t11.*u2p1.*u3p1-M.*t2.*t12.*u2p1.*u3p1-M.*t4.*t12.*u2p1.*u3p1+M.*t2.*t85.*u2p1.*u3p1+t3.*t5.*t25.*t29.*tau1p1+l.*r.*t2.*t3.*t30.*tau3p1+l.*r.*t3.*t4.*t30.*tau3p1+t2.*t3.*t4.*t25.*t29.*tau1p1.*5.0-t3.*t5.*t25.*t28.*t31.*tau2p1+t3.*t5.*t24.*t28.*t30.*tau3p1+I2.*I3.*M.*t4.*u2p1.*u3p1-I2.*I3.*f3p1.*r.*t24.*t31-I3.*M.*f2p1.*r.*t2.*t30-I3.*M.*f2p1.*r.*t4.*t30+I2.*M.*l.*r.*t30.*tau3p1+I3.*M.*t4.*t25.*t29.*tau1p1+I2.*t2.*t3.*t4.*u2p1.*u3p1-I3.*t2.*t3.*t4.*u2p1.*u3p1+M.*t4.*t12.*t25.*u2p1.*u3p1+M.*t4.*t25.*t85.*u2p1.*u3p1-I2.*I3.*M.*t4.*t25.*u2p1.*u3p1.*2.0-I2.*M.*f2p1.*r.*t4.*t25.*t30+I3.*M.*f2p1.*r.*t4.*t25.*t30-I2.*M.*f3p1.*r.*t2.*t24.*t31-I3.*M.*f3p1.*r.*t4.*t24.*t31+I2.*M.*l.*r.*t24.*t28.*tau1p1.*2.0+I3.*M.*l.*r.*t24.*t28.*tau1p1.*2.0-I3.*M.*l.*r.*t24.*t31.*tau2p1-I3.*M.*t4.*t25.*t28.*t31.*tau2p1+I2.*M.*t4.*t24.*t28.*t30.*tau3p1-M.*l.*r.*t30.*t85.*u1p1.*u2p1+I2.*t3.*t5.*t25.*t29.*u2p1.*u3p1-I3.*t3.*t5.*t25.*t29.*u2p1.*u3p1-M.*t4.*t12.*t25.*t29.*u2p1.*u3p1+l.*r.*t2.*t3.*t24.*t28.*tau1p1.*4.0+l.*r.*t3.*t4.*t24.*t28.*tau1p1.*2.0-l.*r.*t2.*t3.*t24.*t31.*tau2p1-l.*r.*t3.*t4.*t24.*t31.*tau2p1-t2.*t3.*t4.*t25.*t28.*t31.*tau2p1.*3.0+t2.*t3.*t4.*t24.*t28.*t30.*tau3p1.*3.0+I1.*I2.*M.*l.*r.*t30.*u1p1.*u2p1+I2.*I3.*M.*l.*r.*t30.*u1p1.*u2p1+I2.*I3.*M.*t4.*t25.*t29.*u2p1.*u3p1+I2.*M.*f1p1.*l.*t4.*t24.*t30.*t31-I3.*M.*f1p1.*l.*t4.*t24.*t30.*t31-I3.*M.*f2p1.*l.*t4.*t24.*t28.*t30.*2.0-I2.*M.*f3p1.*l.*t4.*t25.*t28.*t31.*2.0+I2.*M.*f2p1.*r.*t4.*t25.*t29.*t30-I3.*M.*f2p1.*r.*t4.*t25.*t29.*t30-I2.*M.*f3p1.*r.*t4.*t24.*t25.*t31+I3.*M.*f3p1.*r.*t4.*t24.*t25.*t31+I2.*g.*r.*t2.*t3.*t24.*t30.*t31-I3.*g.*r.*t2.*t3.*t24.*t30.*t31+I1.*l.*r.*t2.*t3.*t30.*u1p1.*u2p1-I2.*l.*r.*t2.*t3.*t30.*u1p1.*u2p1+I1.*l.*r.*t3.*t4.*t30.*u1p1.*u2p1+I3.*l.*r.*t2.*t3.*t30.*u1p1.*u2p1-I2.*l.*r.*t3.*t4.*t30.*u1p1.*u2p1+I3.*l.*r.*t3.*t4.*t30.*u1p1.*u2p1-M.*l.*r.*t12.*t24.*t31.*u1p1.*u3p1-M.*l.*r.*t12.*t24.*t28.*u2p1.*u3p1.*2.0+M.*l.*r.*t24.*t28.*t85.*u2p1.*u3p1.*2.0+I2.*t2.*t3.*t4.*t24.*t30.*t31.*t32+I2.*t2.*t3.*t4.*t24.*t30.*t31.*t33-I3.*t2.*t3.*t4.*t24.*t30.*t31.*t32-I3.*t2.*t3.*t4.*t24.*t30.*t31.*t33+I2.*t2.*t3.*t4.*t25.*t29.*u2p1.*u3p1.*5.0-I3.*t2.*t3.*t4.*t25.*t29.*u2p1.*u3p1.*5.0+I1.*t3.*t5.*t24.*t28.*t30.*u1p1.*u2p1-I2.*t3.*t5.*t24.*t28.*t30.*u1p1.*u2p1+I1.*t3.*t5.*t25.*t28.*t31.*u1p1.*u3p1-I3.*t3.*t5.*t25.*t28.*t31.*u1p1.*u3p1-M.*t4.*t12.*t25.*t28.*t31.*u1p1.*u3p1-M.*t4.*t24.*t28.*t30.*t85.*u1p1.*u2p1+l.*r.*t3.*t4.*t25.*t29.*t30.*tau3p1.*2.0+I2.*M.*f1p1.*r.*t4.*t25.*t28.*t30.*t31-I3.*M.*f1p1.*r.*t4.*t25.*t28.*t30.*t31+I2.*g.*l.*t3.*t4.*t25.*t28.*t30.*t31-I3.*g.*l.*t3.*t4.*t25.*t28.*t30.*t31+I2.*l.*r.*t3.*t4.*t25.*t30.*u1p1.*u2p1-I3.*l.*r.*t3.*t4.*t25.*t30.*u1p1.*u2p1+I1.*l.*r.*t2.*t3.*t24.*t31.*u1p1.*u3p1+I2.*l.*r.*t2.*t3.*t24.*t31.*u1p1.*u3p1+I1.*l.*r.*t3.*t4.*t24.*t31.*u1p1.*u3p1-I3.*l.*r.*t2.*t3.*t24.*t31.*u1p1.*u3p1+I2.*l.*r.*t2.*t3.*t24.*t28.*u2p1.*u3p1.*4.0-I3.*l.*r.*t2.*t3.*t24.*t28.*u2p1.*u3p1.*4.0+I2.*l.*r.*t3.*t4.*t24.*t28.*u2p1.*u3p1.*2.0-I3.*l.*r.*t3.*t4.*t24.*t28.*u2p1.*u3p1.*2.0+I1.*t2.*t3.*t4.*t24.*t28.*t30.*u1p1.*u2p1.*3.0-I2.*t2.*t3.*t4.*t24.*t28.*t30.*u1p1.*u2p1.*3.0+I3.*t2.*t3.*t4.*t24.*t28.*t30.*u1p1.*u2p1.*2.0+I1.*t2.*t3.*t4.*t25.*t28.*t31.*u1p1.*u3p1.*3.0+I2.*t2.*t3.*t4.*t25.*t28.*t31.*u1p1.*u3p1.*2.0-I3.*t2.*t3.*t4.*t25.*t28.*t31.*u1p1.*u3p1.*3.0+l.*r.*t3.*t4.*t24.*t25.*t28.*t29.*tau1p1.*2.0-l.*r.*t3.*t4.*t24.*t25.*t29.*t31.*tau2p1.*2.0+I1.*I3.*M.*l.*r.*t24.*t31.*u1p1.*u3p1+I2.*I3.*M.*l.*r.*t24.*t31.*u1p1.*u3p1+I1.*I2.*M.*t4.*t24.*t28.*t30.*u1p1.*u2p1+I1.*I3.*M.*t4.*t25.*t28.*t31.*u1p1.*u3p1+I2.*l.*r.*t3.*t4.*t25.*t28.*t30.*t31.*t32+I2.*l.*r.*t3.*t4.*t25.*t28.*t30.*t31.*t33-I3.*l.*r.*t3.*t4.*t25.*t28.*t30.*t31.*t32-I3.*l.*r.*t3.*t4.*t25.*t28.*t30.*t31.*t33+I1.*l.*r.*t3.*t4.*t25.*t29.*t30.*u1p1.*u2p1.*2.0-I2.*l.*r.*t3.*t4.*t25.*t29.*t30.*u1p1.*u2p1.*3.0+I3.*l.*r.*t3.*t4.*t25.*t29.*t30.*u1p1.*u2p1+I2.*l.*r.*t3.*t4.*t24.*t25.*t31.*u1p1.*u3p1-I3.*l.*r.*t3.*t4.*t24.*t25.*t31.*u1p1.*u3p1+I1.*l.*r.*t3.*t4.*t24.*t25.*t29.*t31.*u1p1.*u3p1.*2.0-I3.*l.*r.*t3.*t4.*t24.*t25.*t29.*t31.*u1p1.*u3p1.*2.0+I2.*l.*r.*t3.*t4.*t24.*t25.*t28.*t29.*u2p1.*u3p1.*2.0-I3.*l.*r.*t3.*t4.*t24.*t25.*t28.*t29.*u2p1.*u3p1.*2.0)).*(1.0./8.0)).*(dt.*t28.*t37.*t48.*t57.*(-1.0./8.0)+dt.*t31.*t48.*t66.*t101.*(1.0./8.0)-dt.*t28.*t30.*t37.*t58.*t66.*(1.0./8.0)+dt.*t28.*t37.*t58.*t66.*t69.*(1.0./8.0)+dt.*t28.*t30.*t37.*t48.*t57.*t69.*(1.0./8.0))+(dt.*t100.*(-M.*l.*r.*t30.*t85.*u2p1+I1.*I2.*M.*l.*r.*t30.*u2p1+I2.*I3.*M.*l.*r.*t30.*u2p1+I1.*l.*r.*t2.*t3.*t30.*u2p1-I2.*l.*r.*t2.*t3.*t30.*u2p1+I1.*l.*r.*t3.*t4.*t30.*u2p1+I3.*l.*r.*t2.*t3.*t30.*u2p1-I2.*l.*r.*t3.*t4.*t30.*u2p1+I3.*l.*r.*t3.*t4.*t30.*u2p1-M.*l.*r.*t12.*t24.*t31.*u3p1+I1.*t3.*t5.*t24.*t28.*t30.*u2p1-I2.*t3.*t5.*t24.*t28.*t30.*u2p1+I1.*t3.*t5.*t25.*t28.*t31.*u3p1-I3.*t3.*t5.*t25.*t28.*t31.*u3p1-M.*t4.*t12.*t25.*t28.*t31.*u3p1-M.*t4.*t24.*t28.*t30.*t85.*u2p1+I1.*I3.*M.*l.*r.*t24.*t31.*u3p1+I2.*I3.*M.*l.*r.*t24.*t31.*u3p1+I1.*I2.*M.*t4.*t24.*t28.*t30.*u2p1+I1.*I3.*M.*t4.*t25.*t28.*t31.*u3p1+I2.*l.*r.*t3.*t4.*t25.*t30.*u2p1-I3.*l.*r.*t3.*t4.*t25.*t30.*u2p1+I1.*l.*r.*t2.*t3.*t24.*t31.*u3p1+I2.*l.*r.*t2.*t3.*t24.*t31.*u3p1+I1.*l.*r.*t3.*t4.*t24.*t31.*u3p1-I3.*l.*r.*t2.*t3.*t24.*t31.*u3p1+I1.*t2.*t3.*t4.*t24.*t28.*t30.*u2p1.*3.0-I2.*t2.*t3.*t4.*t24.*t28.*t30.*u2p1.*3.0+I3.*t2.*t3.*t4.*t24.*t28.*t30.*u2p1.*2.0+I1.*t2.*t3.*t4.*t25.*t28.*t31.*u3p1.*3.0+I2.*t2.*t3.*t4.*t25.*t28.*t31.*u3p1.*2.0-I3.*t2.*t3.*t4.*t25.*t28.*t31.*u3p1.*3.0+I1.*l.*r.*t3.*t4.*t25.*t29.*t30.*u2p1.*2.0-I2.*l.*r.*t3.*t4.*t25.*t29.*t30.*u2p1.*3.0+I3.*l.*r.*t3.*t4.*t25.*t29.*t30.*u2p1+I2.*l.*r.*t3.*t4.*t24.*t25.*t31.*u3p1-I3.*l.*r.*t3.*t4.*t24.*t25.*t31.*u3p1+I1.*l.*r.*t3.*t4.*t24.*t25.*t29.*t31.*u3p1.*2.0-I3.*l.*r.*t3.*t4.*t24.*t25.*t29.*t31.*u3p1.*2.0).*(1.0./8.0)-1.0./2.0).*(t57.*t58+t48.*t66.*t69)-(u2.*(1.0./2.0)+u2p1.*(1.0./2.0)+dt.*(t84.*(I1.*I3.*tau2+I1.*t12.*u1.*u3-I3.*t10.*u1.*u3+t3.*t5.*t7.*tau2-I1.*I3.*f3.*l+I1.*M.*t2.*tau2+I3.*M.*t4.*tau2-I1.*M.*f3.*l.*t2-I3.*M.*f3.*l.*t4+I1.*I3.*f1.*r.*t13+I1.*M.*t4.*t7.*tau2+I1.*t3.*t11.*u1.*u3-M.*t2.*t10.*u1.*u3+M.*t4.*t12.*u1.*u3+t2.*t3.*t4.*t7.*tau2-t3.*t5.*t7.*t9.*tau2-t2.*t3.*t4.*t7.*t9.*tau2-t3.*t5.*t7.*t8.*t16.*tau1-t3.*t5.*t6.*t13.*t16.*tau3+I1.*I3.*M.*g.*l.*t13+I1.*I3.*M.*t2.*u1.*u3.*2.0-I1.*I3.*M.*t4.*u1.*u3-I1.*M.*f3.*l.*t4.*t7+I3.*M.*f3.*l.*t4.*t7-I1.*I3.*f3.*r.*t6.*t8+I1.*M.*f1.*r.*t2.*t13+I3.*M.*f1.*r.*t4.*t13-I3.*M.*t4.*t7.*t9.*tau2+I1.*g.*l.*t2.*t3.*t13+I3.*g.*l.*t3.*t4.*t13+I3.*t2.*t3.*t4.*u1.*u3-I1.*t3.*t5.*t7.*u1.*u3+I3.*t3.*t5.*t7.*u1.*u3-M.*t4.*t7.*t10.*u1.*u3+I1.*I3.*M.*l.*r.*t13.*t14+I1.*I3.*M.*l.*r.*t13.*t15+I1.*I3.*M.*t4.*t7.*u1.*u3-I1.*M.*f3.*l.*t4.*t7.*t9.*2.0-I1.*M.*f3.*r.*t2.*t6.*t8.*3.0-I3.*M.*f3.*r.*t4.*t6.*t8+I1.*M.*l.*r.*t6.*t8.*tau2.*2.0-I3.*M.*l.*r.*t6.*t16.*tau1-I3.*M.*t4.*t7.*t8.*t16.*tau1-I1.*M.*t4.*t6.*t13.*t16.*tau3+I1.*l.*r.*t2.*t3.*t13.*t14+I1.*l.*r.*t2.*t3.*t13.*t15+I3.*l.*r.*t3.*t4.*t13.*t14+I3.*l.*r.*t3.*t4.*t13.*t15+I1.*t3.*t5.*t7.*t9.*u1.*u3-I3.*t3.*t5.*t7.*t9.*u1.*u3-M.*t4.*t7.*t9.*t12.*u1.*u3-l.*r.*t2.*t3.*t6.*t16.*tau1-l.*r.*t3.*t4.*t6.*t16.*tau1-t2.*t3.*t4.*t7.*t8.*t16.*tau1.*3.0-t2.*t3.*t4.*t6.*t13.*t16.*tau3+I1.*I3.*M.*t4.*t7.*t9.*u1.*u3+I1.*M.*f1.*l.*t4.*t6.*t8.*t13.*2.0-I1.*M.*f2.*l.*t4.*t6.*t13.*t16+I3.*M.*f2.*l.*t4.*t6.*t13.*t16-I1.*M.*f3.*r.*t4.*t6.*t7.*t8+I3.*M.*f3.*r.*t4.*t6.*t7.*t8+I1.*M.*f1.*r.*t4.*t7.*t9.*t13-I3.*M.*f1.*r.*t4.*t7.*t9.*t13+I1.*g.*l.*t3.*t4.*t7.*t9.*t13-I3.*g.*l.*t3.*t4.*t7.*t9.*t13+I1.*g.*r.*t2.*t3.*t6.*t8.*t13.*2.0-M.*l.*r.*t6.*t8.*t10.*u1.*u3.*2.0+M.*l.*r.*t6.*t12.*t16.*u2.*u3+I1.*t2.*t3.*t4.*t6.*t8.*t13.*t14.*2.0+I1.*t2.*t3.*t4.*t6.*t8.*t13.*t15.*2.0+I1.*t2.*t3.*t4.*t7.*t9.*u1.*u3.*3.0-I3.*t2.*t3.*t4.*t7.*t9.*u1.*u3-I2.*t3.*t5.*t7.*t8.*t16.*u2.*u3-I1.*t3.*t5.*t6.*t13.*t16.*u1.*u2+I3.*t3.*t5.*t7.*t8.*t16.*u2.*u3+I2.*t3.*t5.*t6.*t13.*t16.*u1.*u2-M.*t4.*t6.*t10.*t13.*t16.*u1.*u2+M.*t4.*t7.*t8.*t12.*t16.*u2.*u3+l.*r.*t3.*t4.*t6.*t7.*t8.*tau2.*2.0-I1.*M.*f2.*r.*t4.*t7.*t8.*t13.*t16+I3.*M.*f2.*r.*t4.*t7.*t8.*t13.*t16+I1.*l.*r.*t3.*t4.*t7.*t9.*t13.*t14+I1.*l.*r.*t3.*t4.*t7.*t9.*t13.*t15-I3.*l.*r.*t3.*t4.*t7.*t9.*t13.*t14-I3.*l.*r.*t3.*t4.*t7.*t9.*t13.*t15+I1.*l.*r.*t2.*t3.*t6.*t8.*u1.*u3.*3.0+I3.*l.*r.*t3.*t4.*t6.*t8.*u1.*u3-I2.*l.*r.*t2.*t3.*t6.*t16.*u2.*u3+I3.*l.*r.*t2.*t3.*t6.*t16.*u2.*u3-I2.*l.*r.*t3.*t4.*t6.*t16.*u2.*u3+I3.*l.*r.*t3.*t4.*t6.*t16.*u2.*u3-I2.*t2.*t3.*t4.*t7.*t8.*t16.*u2.*u3.*3.0+I3.*t2.*t3.*t4.*t7.*t8.*t16.*u2.*u3.*3.0+I2.*t2.*t3.*t4.*t6.*t13.*t16.*u1.*u2-I3.*t2.*t3.*t4.*t6.*t13.*t16.*u1.*u2-l.*r.*t3.*t4.*t6.*t7.*t8.*t9.*tau2.*2.0-l.*r.*t3.*t4.*t6.*t7.*t9.*t16.*tau1.*2.0-l.*r.*t3.*t4.*t7.*t8.*t13.*t16.*tau3.*2.0+I1.*I3.*M.*l.*r.*t6.*t8.*u1.*u3.*3.0-I2.*I3.*M.*l.*r.*t6.*t16.*u2.*u3+I1.*I2.*M.*t4.*t6.*t13.*t16.*u1.*u2-I2.*I3.*M.*t4.*t7.*t8.*t16.*u2.*u3-I1.*l.*r.*t3.*t4.*t6.*t7.*t8.*u1.*u3+I3.*l.*r.*t3.*t4.*t6.*t7.*t8.*u1.*u3+I1.*l.*r.*t3.*t4.*t6.*t7.*t8.*t9.*u1.*u3.*2.0-I3.*l.*r.*t3.*t4.*t6.*t7.*t8.*t9.*u1.*u3.*2.0-I2.*l.*r.*t3.*t4.*t6.*t7.*t9.*t16.*u2.*u3.*2.0+I3.*l.*r.*t3.*t4.*t6.*t7.*t9.*t16.*u2.*u3.*2.0-I1.*l.*r.*t3.*t4.*t7.*t8.*t13.*t16.*u1.*u2+I2.*l.*r.*t3.*t4.*t7.*t8.*t13.*t16.*u1.*u2.*2.0-I3.*l.*r.*t3.*t4.*t7.*t8.*t13.*t16.*u1.*u2)-t100.*(I1.*I3.*tau2p1+I1.*t12.*u1p1.*u3p1-I3.*t10.*u1p1.*u3p1+t3.*t5.*t25.*tau2p1-I1.*I3.*f3p1.*l+I1.*M.*t2.*tau2p1+I3.*M.*t4.*tau2p1-I1.*M.*f3p1.*l.*t2-I3.*M.*f3p1.*l.*t4+I1.*I3.*f1p1.*r.*t30+I1.*M.*t4.*t25.*tau2p1+I1.*t3.*t11.*u1p1.*u3p1-M.*t2.*t10.*u1p1.*u3p1+M.*t4.*t12.*u1p1.*u3p1+t2.*t3.*t4.*t25.*tau2p1-t3.*t5.*t25.*t29.*tau2p1-t2.*t3.*t4.*t25.*t29.*tau2p1-t3.*t5.*t25.*t28.*t31.*tau1p1-t3.*t5.*t24.*t30.*t31.*tau3p1+I1.*I3.*M.*g.*l.*t30+I1.*I3.*M.*t2.*u1p1.*u3p1.*2.0-I1.*I3.*M.*t4.*u1p1.*u3p1-I1.*M.*f3p1.*l.*t4.*t25+I3.*M.*f3p1.*l.*t4.*t25-I1.*I3.*f3p1.*r.*t24.*t28+I1.*M.*f1p1.*r.*t2.*t30+I3.*M.*f1p1.*r.*t4.*t30-I3.*M.*t4.*t25.*t29.*tau2p1+I1.*g.*l.*t2.*t3.*t30+I3.*g.*l.*t3.*t4.*t30+I3.*t2.*t3.*t4.*u1p1.*u3p1-I1.*t3.*t5.*t25.*u1p1.*u3p1+I3.*t3.*t5.*t25.*u1p1.*u3p1-M.*t4.*t10.*t25.*u1p1.*u3p1+I1.*I3.*M.*l.*r.*t30.*t32+I1.*I3.*M.*l.*r.*t30.*t33+I1.*I3.*M.*t4.*t25.*u1p1.*u3p1-I1.*M.*f3p1.*l.*t4.*t25.*t29.*2.0-I1.*M.*f3p1.*r.*t2.*t24.*t28.*3.0-I3.*M.*f3p1.*r.*t4.*t24.*t28-I3.*M.*l.*r.*t24.*t31.*tau1p1+I1.*M.*l.*r.*t24.*t28.*tau2p1.*2.0-I3.*M.*t4.*t25.*t28.*t31.*tau1p1-I1.*M.*t4.*t24.*t30.*t31.*tau3p1+I1.*l.*r.*t2.*t3.*t30.*t32+I1.*l.*r.*t2.*t3.*t30.*t33+I3.*l.*r.*t3.*t4.*t30.*t32+I3.*l.*r.*t3.*t4.*t30.*t33+I1.*t3.*t5.*t25.*t29.*u1p1.*u3p1-I3.*t3.*t5.*t25.*t29.*u1p1.*u3p1-M.*t4.*t12.*t25.*t29.*u1p1.*u3p1-l.*r.*t2.*t3.*t24.*t31.*tau1p1-l.*r.*t3.*t4.*t24.*t31.*tau1p1-t2.*t3.*t4.*t25.*t28.*t31.*tau1p1.*3.0-t2.*t3.*t4.*t24.*t30.*t31.*tau3p1+I1.*I3.*M.*t4.*t25.*t29.*u1p1.*u3p1+I1.*M.*f1p1.*l.*t4.*t24.*t28.*t30.*2.0-I1.*M.*f2p1.*l.*t4.*t24.*t30.*t31+I3.*M.*f2p1.*l.*t4.*t24.*t30.*t31+I1.*M.*f1p1.*r.*t4.*t25.*t29.*t30-I3.*M.*f1p1.*r.*t4.*t25.*t29.*t30-I1.*M.*f3p1.*r.*t4.*t24.*t25.*t28+I3.*M.*f3p1.*r.*t4.*t24.*t25.*t28+I1.*g.*l.*t3.*t4.*t25.*t29.*t30-I3.*g.*l.*t3.*t4.*t25.*t29.*t30+I1.*g.*r.*t2.*t3.*t24.*t28.*t30.*2.0-M.*l.*r.*t10.*t24.*t28.*u1p1.*u3p1.*2.0+M.*l.*r.*t12.*t24.*t31.*u2p1.*u3p1+I1.*t2.*t3.*t4.*t24.*t28.*t30.*t32.*2.0+I1.*t2.*t3.*t4.*t24.*t28.*t30.*t33.*2.0+I1.*t2.*t3.*t4.*t25.*t29.*u1p1.*u3p1.*3.0-I3.*t2.*t3.*t4.*t25.*t29.*u1p1.*u3p1-I1.*t3.*t5.*t24.*t30.*t31.*u1p1.*u2p1+I2.*t3.*t5.*t24.*t30.*t31.*u1p1.*u2p1-I2.*t3.*t5.*t25.*t28.*t31.*u2p1.*u3p1+I3.*t3.*t5.*t25.*t28.*t31.*u2p1.*u3p1-M.*t4.*t10.*t24.*t30.*t31.*u1p1.*u2p1+M.*t4.*t12.*t25.*t28.*t31.*u2p1.*u3p1+l.*r.*t3.*t4.*t24.*t25.*t28.*tau2p1.*2.0-I1.*M.*f2p1.*r.*t4.*t25.*t28.*t30.*t31+I3.*M.*f2p1.*r.*t4.*t25.*t28.*t30.*t31+I1.*l.*r.*t3.*t4.*t25.*t29.*t30.*t32+I1.*l.*r.*t3.*t4.*t25.*t29.*t30.*t33-I3.*l.*r.*t3.*t4.*t25.*t29.*t30.*t32-I3.*l.*r.*t3.*t4.*t25.*t29.*t30.*t33+I1.*l.*r.*t2.*t3.*t24.*t28.*u1p1.*u3p1.*3.0+I3.*l.*r.*t3.*t4.*t24.*t28.*u1p1.*u3p1-I2.*l.*r.*t2.*t3.*t24.*t31.*u2p1.*u3p1+I3.*l.*r.*t2.*t3.*t24.*t31.*u2p1.*u3p1-I2.*l.*r.*t3.*t4.*t24.*t31.*u2p1.*u3p1+I3.*l.*r.*t3.*t4.*t24.*t31.*u2p1.*u3p1+I2.*t2.*t3.*t4.*t24.*t30.*t31.*u1p1.*u2p1-I3.*t2.*t3.*t4.*t24.*t30.*t31.*u1p1.*u2p1-I2.*t2.*t3.*t4.*t25.*t28.*t31.*u2p1.*u3p1.*3.0+I3.*t2.*t3.*t4.*t25.*t28.*t31.*u2p1.*u3p1.*3.0-l.*r.*t3.*t4.*t24.*t25.*t29.*t31.*tau1p1.*2.0-l.*r.*t3.*t4.*t24.*t25.*t28.*t29.*tau2p1.*2.0-l.*r.*t3.*t4.*t25.*t28.*t30.*t31.*tau3p1.*2.0+I1.*I3.*M.*l.*r.*t24.*t28.*u1p1.*u3p1.*3.0-I2.*I3.*M.*l.*r.*t24.*t31.*u2p1.*u3p1+I1.*I2.*M.*t4.*t24.*t30.*t31.*u1p1.*u2p1-I2.*I3.*M.*t4.*t25.*t28.*t31.*u2p1.*u3p1-I1.*l.*r.*t3.*t4.*t24.*t25.*t28.*u1p1.*u3p1+I3.*l.*r.*t3.*t4.*t24.*t25.*t28.*u1p1.*u3p1-I1.*l.*r.*t3.*t4.*t25.*t28.*t30.*t31.*u1p1.*u2p1+I2.*l.*r.*t3.*t4.*t25.*t28.*t30.*t31.*u1p1.*u2p1.*2.0+I1.*l.*r.*t3.*t4.*t24.*t25.*t28.*t29.*u1p1.*u3p1.*2.0-I3.*l.*r.*t3.*t4.*t25.*t28.*t30.*t31.*u1p1.*u2p1-I3.*l.*r.*t3.*t4.*t24.*t25.*t28.*t29.*u1p1.*u3p1.*2.0-I2.*l.*r.*t3.*t4.*t24.*t25.*t29.*t31.*u2p1.*u3p1.*2.0+I3.*l.*r.*t3.*t4.*t24.*t25.*t29.*t31.*u2p1.*u3p1.*2.0)).*(1.0./8.0)).*(dt.*t28.*t37.*t48.*t66.*(1.0./8.0)+dt.*t31.*t48.*t57.*t101.*(1.0./8.0)-dt.*t28.*t30.*t37.*t57.*t58.*(1.0./8.0)+dt.*t28.*t37.*t57.*t58.*t69.*(1.0./8.0)-dt.*t28.*t30.*t37.*t48.*t66.*t69.*(1.0./8.0))+dt.*t100.*(t58.*t66-t48.*t57.*t69).*(I1.*t12.*u3p1-I3.*t10.*u3p1+I1.*t3.*t11.*u3p1-M.*t2.*t10.*u3p1+M.*t4.*t12.*u3p1+I1.*I3.*M.*t2.*u3p1.*2.0-I1.*I3.*M.*t4.*u3p1+I3.*t2.*t3.*t4.*u3p1-I1.*t3.*t5.*t25.*u3p1+I3.*t3.*t5.*t25.*u3p1-M.*t4.*t10.*t25.*u3p1+I1.*I3.*M.*t4.*t25.*u3p1+I1.*t3.*t5.*t25.*t29.*u3p1-I3.*t3.*t5.*t25.*t29.*u3p1-M.*t4.*t12.*t25.*t29.*u3p1+I1.*I3.*M.*t4.*t25.*t29.*u3p1-M.*l.*r.*t10.*t24.*t28.*u3p1.*2.0+I1.*t2.*t3.*t4.*t25.*t29.*u3p1.*3.0-I3.*t2.*t3.*t4.*t25.*t29.*u3p1-I1.*t3.*t5.*t24.*t30.*t31.*u2p1+I2.*t3.*t5.*t24.*t30.*t31.*u2p1-M.*t4.*t10.*t24.*t30.*t31.*u2p1+I1.*I3.*M.*l.*r.*t24.*t28.*u3p1.*3.0+I1.*I2.*M.*t4.*t24.*t30.*t31.*u2p1+I1.*l.*r.*t2.*t3.*t24.*t28.*u3p1.*3.0+I3.*l.*r.*t3.*t4.*t24.*t28.*u3p1+I2.*t2.*t3.*t4.*t24.*t30.*t31.*u2p1-I3.*t2.*t3.*t4.*t24.*t30.*t31.*u2p1-I1.*l.*r.*t3.*t4.*t24.*t25.*t28.*u3p1+I3.*l.*r.*t3.*t4.*t24.*t25.*t28.*u3p1-I1.*l.*r.*t3.*t4.*t25.*t28.*t30.*t31.*u2p1+I2.*l.*r.*t3.*t4.*t25.*t28.*t30.*t31.*u2p1.*2.0+I1.*l.*r.*t3.*t4.*t24.*t25.*t28.*t29.*u3p1.*2.0-I3.*l.*r.*t3.*t4.*t25.*t28.*t30.*t31.*u2p1-I3.*l.*r.*t3.*t4.*t24.*t25.*t28.*t29.*u3p1.*2.0).*(1.0./8.0)+dt.*t31.*t48.*t69.*t271.*(1.0./8.0)+dt.*t48.*t100.*t101.*(-I2.*t10.*u2p1+I1.*t85.*u2p1-I1.*t3.*t5.*u2p1+I2.*t3.*t5.*u2p1+I1.*t3.*t11.*u2p1-M.*t2.*t10.*u2p1-M.*t4.*t10.*u2p1+M.*t4.*t85.*u2p1+I1.*I2.*M.*t2.*u2p1.*2.0+I2.*t2.*t3.*t4.*u2p1+I1.*t3.*t5.*t25.*u2p1-I2.*t3.*t5.*t25.*u2p1+M.*t4.*t10.*t25.*u2p1-I1.*I2.*M.*t4.*t25.*u2p1-M.*t4.*t10.*t25.*t29.*u2p1-M.*t4.*t25.*t29.*t85.*u2p1+I1.*I2.*M.*t4.*t25.*t29.*u2p1.*2.0-M.*l.*r.*t10.*t24.*t28.*u2p1.*2.0+I1.*t2.*t3.*t4.*t25.*t29.*u2p1.*3.0-I2.*t2.*t3.*t4.*t25.*t29.*u2p1-I1.*t3.*t5.*t24.*t30.*t31.*u3p1+I3.*t3.*t5.*t24.*t30.*t31.*u3p1-M.*t4.*t10.*t24.*t30.*t31.*u3p1+I1.*I2.*M.*l.*r.*t24.*t28.*u2p1.*3.0+I1.*I3.*M.*t4.*t24.*t30.*t31.*u3p1+I1.*l.*r.*t2.*t3.*t24.*t28.*u2p1.*3.0-I1.*l.*r.*t3.*t4.*t24.*t28.*u2p1+I2.*l.*r.*t3.*t4.*t24.*t28.*u2p1.*2.0-I2.*t2.*t3.*t4.*t24.*t30.*t31.*u3p1+I3.*t2.*t3.*t4.*t24.*t30.*t31.*u3p1+I1.*l.*r.*t3.*t4.*t24.*t25.*t28.*u2p1-I2.*l.*r.*t3.*t4.*t24.*t25.*t28.*u2p1+I1.*l.*r.*t3.*t4.*t24.*t25.*t28.*t29.*u2p1-I2.*l.*r.*t3.*t4.*t24.*t25.*t28.*t29.*u2p1-I1.*l.*r.*t3.*t4.*t25.*t28.*t30.*t31.*u3p1-I2.*l.*r.*t3.*t4.*t25.*t28.*t30.*t31.*u3p1+I3.*l.*r.*t3.*t4.*t25.*t28.*t30.*t31.*u3p1.*2.0).*(1.0./8.0)-dt.*t28.*t37.*t58.*t101.*t271.*(1.0./8.0));
