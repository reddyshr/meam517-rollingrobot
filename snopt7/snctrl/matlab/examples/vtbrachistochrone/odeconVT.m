function [F,J]=odeconVT (snStat,curPhs,nPhs,nY,nU,nP,nNodes,...
			 dvar,pvar,needF,needJ)

F = 0;
J = 0;

if (needF > 0)
  for jt = 1:nNodes,
    F(1,jt) = dvar(3,jt)*cos(dvar(1+nY,jt))*pvar(1);
    F(2,jt) = dvar(3,jt)*sin(dvar(1+nY,jt))*pvar(1);
    F(3,jt) = sin(dvar(1+nY,jt))*pvar(1);
    F(4,jt) = pvar(1);
  end
end


if (needJ > 0)
  for jt = 1:nNodes,
    J(1,jt) = 0;
    J(5,jt) = 0;
    J(9,jt) = cos(dvar(1+nY,jt))*pvar(1);
    J(13,jt) = 0;
    J(17,jt) = -dvar(3,jt)*sin(dvar(1+nY,jt))*pvar(1);
    J(21,jt) = dvar(3,jt)*cos(dvar(1+nY,jt));

    J(2,jt) = 0;
    J(6,jt) = 0;
    J(10,jt) = sin(dvar(1+nY,jt))*pvar(1);
    J(14,jt) = 0;
    J(18,jt) = dvar(3,jt)*cos(dvar(1+nY,jt))*pvar(1);
    J(22,jt) = dvar(3,jt)*sin(dvar(1+nY,jt));

    J(3,jt) = 0;
    J(7,jt) = 0;
    J(11,jt) = 0;
    J(15,jt) = 0;
    J(19,jt) = cos(dvar(1+nY,jt))*pvar(1);
    J(23,jt) = sin(dvar(1+nY,jt));

    J(4,jt) = 0;
    J(8,jt) = 0;
    J(12,jt) = 0;
    J(16,jt) = 0;
    J(20,jt) = 0;
    J(24,jt) = 1;
  end
end


