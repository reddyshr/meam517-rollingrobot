function control_des = getContinuousInputTrajectory(u_col, N, dt, timestep)
    
    %% Piecewise Linear Input Spline
    
    numCols = length(u_col);
    numRows = length(u_col(:,1));
    
    control_des = zeros(numRows, ((N)*ceil(dt/timestep)));
    
    tdiff = 0;
    seg_cnt = 1;
    
    for i = 1:length(control_des)
        if (tdiff >= dt) 
            seg_cnt = seg_cnt + 1;
            tdiff = 0;
        end
        if (seg_cnt >= numCols)
            break;
        end        
        vect = (tdiff/dt) * (u_col(:,seg_cnt+1) - u_col(:,seg_cnt)) +  u_col(:,seg_cnt);
        
        control_des(:,i) = vect;
        
        tdiff = tdiff + timestep;
        
    end
    
end