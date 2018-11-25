function [F,J]=odeconCM (snStat,curPhs,nPhs,nY,nU,nP,nNodes,...
			 dvar,pvar,needF,needJ)

F = 0;
J = 0;

if (needF > 0)
  for jt = 1:nNodes,
    F(1,jt) = dvar(1+nY,jt)*(10*dvar(2,jt) - dvar(1,jt));
    F(2,jt) = dvar(1+nY,jt)*(dvar(1,jt) - 10*dvar(2,jt)) - (1-dvar(1+nY,jt))*dvar(2,jt);
  end
end


if (needJ > 0)
  for jt = 1:nNodes,
    J(1,jt) = -dvar(1+nY,jt);
    J(3,jt) = 10*dvar(1+nY,jt);
    J(5,jt) = 10*dvar(2,jt) - dvar(1,jt);
    J(7,jt) = 0;

    J(2,jt) = dvar(1+nY,jt);
    J(4,jt) = -10*dvar(1+nY,jt) - (1-dvar(1+nY,jt));
    J(6,jt) = dvar(1,jt) - 10*dvar(2,jt) + dvar(2,jt);
    J(8,jt) = 0;
  end
end


