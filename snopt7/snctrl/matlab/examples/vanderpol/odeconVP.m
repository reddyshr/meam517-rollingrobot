function [F,J]=odeconVP (snStat,curPhs,nPhs,nY,nU,nP,nNodes,...
			 dvar,pvar,needF,needJ)

F = 0;
J = 0;

if (needF > 0)
  for jt = 1:nNodes,
    F(1,jt) = (dvar(2,jt)^2 + dvar(3,jt)^2 + dvar(1+nY,jt)^2)/2d+0;
    F(2,jt) = -dvar(3,jt) + (1 - dvar(3,jt)^2)*dvar(2,jt) + dvar(1+nY,jt);
    F(3,jt) = dvar(2,jt);
  end
end


if (needJ > 0)
  for jt = 1:nNodes,
    J(1,jt) = 0;
    J(4,jt) = dvar(2,jt);
    J(7,jt) = dvar(3,jt);
    J(10,jt) = dvar(1+nY,jt);

    J(2,jt) = 0;
    J(5,jt) = 1 - dvar(3,jt)^2;
    J(8,jt) = -1 - 2*dvar(2,jt)*dvar(3,jt);
    J(11,jt) = 1;

    J(3,jt) = 0;
    J(6,jt) = 1;
    J(9,jt) = 0;
    J(12,jt) = 0;
  end
end


