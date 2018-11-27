function state_des = getContinuousStateTrajectory(x_col, u_col, N, dt, timestep,M,l,r,g,I1,I2,I3)    
    %% Piecewise Linear Input Spline
    
    numCols = length(x_col);
    numRows = length(x_col(:,1));
    
    state_des = zeros(numRows, ((N)*ceil(dt/timestep)));
    
    tdiff = 0;
    seg_cnt = 1;
    firstTime = true;
    for i = 1:length(state_des)
        if (firstTime)
            [f_k, f_kp1, f_kphalf] = getSplineParams(x_col, u_col, dt, seg_cnt, M,l,r,g,I1,I2,I3);
            firstTime = false;
        elseif (tdiff >= dt) 
            seg_cnt = seg_cnt + 1;
            tdiff = 0;
            if (seg_cnt >= numCols)
                break;
            end    
            [f_k, f_kp1, f_kphalf] = getSplineParams(x_col, u_col, dt, seg_cnt, M,l,r,g,I1,I2,I3);
        end
          
        if (seg_cnt == 29)
            hiya = 1;
        end
        
        vect = x_col(:,seg_cnt) + f_k*(tdiff/dt) +  ... 
               (1/2)*(-3*f_k + 4*f_kphalf - f_kp1)*((tdiff/dt)^2) + ... 
               (1/3)*(2*f_k - 4*f_kphalf + 2*f_kp1)*((tdiff/dt)^3);
        
        state_des(:,i) = vect;
        
        tdiff = tdiff + timestep;
        
    end
    
end

function [f_k, f_kp1, f_kphalf] = getSplineParams(x_col, u_col, dt, seg_cnt, M,l,r,g,I1,I2,I3)

 f_k = F(x_col(1,seg_cnt), x_col(2,seg_cnt), x_col(3,seg_cnt), x_col(4,seg_cnt), ... 
                x_col(5,seg_cnt), x_col(6,seg_cnt), x_col(7,seg_cnt), x_col(8,seg_cnt), ...
                 u_col(1,seg_cnt), u_col(2,seg_cnt), u_col(3,seg_cnt), u_col(4,seg_cnt), ... 
                 u_col(5,seg_cnt), u_col(6,seg_cnt), M,l,r,g,I1,I2,I3);
 f_kp1 = F(x_col(1,seg_cnt+1), x_col(2,seg_cnt+1), x_col(3,seg_cnt+1), x_col(4,seg_cnt+1), ... 
        x_col(5,seg_cnt+1), x_col(6,seg_cnt+1), x_col(7,seg_cnt+1), x_col(8,seg_cnt+1), ...
         u_col(1,seg_cnt+1), u_col(2,seg_cnt+1), u_col(3,seg_cnt+1), u_col(4,seg_cnt+1), ... 
         u_col(5,seg_cnt+1), u_col(6,seg_cnt+1), M,l,r,g,I1,I2,I3);

 x_iphalf = 0.5*(x_col(:,seg_cnt) + x_col(:,seg_cnt+1)) - (dt/8)*(f_kp1 - f_k);
 u_iphalf = 0.5*(u_col(:,seg_cnt) + u_col(:,seg_cnt+1));

 f_kphalf = F(x_iphalf(1), x_iphalf(2), x_iphalf(3), x_iphalf(4), ...
              x_iphalf(5), x_iphalf(6), x_iphalf(7), x_iphalf(8), ...
              u_iphalf(1), u_iphalf(2), u_iphalf(3), u_iphalf(4), ...
              u_iphalf(5), u_iphalf(6), M,l,r,g,I1,I2,I3);

end