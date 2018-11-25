function [F,J]=odeconBR (snStat,curPhs,nPhs,nY,nU,nP,nNodes,...
			 dvar,pvar,needF,needJ)

F = 0;
J = 0;

if (needF > 0)
  for jt = 1:nNodes,
    F(1,jt) = sqrt( (1+dvar(1+nY,jt)^2)/(1-dvar(2,jt)) );
    F(2,jt) = dvar(1+nY,jt);
  end
end


if (needJ > 0)
  for jt = 1:nNodes,
    J(1,jt) = 0;
    J(3,jt) = sqrt((1-dvar(2,jt))/(1+dvar(nY+1,jt)^2))*(1+ ...
						  dvar(nY+1,jt)^2)*0.5 /((1-dvar(2,jt))^2);
    J(5,jt) = sqrt((1-dvar(2,jt))/(1+dvar(nY+1,jt)^2))*dvar(nY+1,jt)/(1-dvar(2,jt));

    J(2,jt) = 0;
    J(4,jt) = 0;
    J(6,jt) = 1;
  end
end


