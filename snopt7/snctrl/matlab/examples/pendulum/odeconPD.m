function [F,J]=odeconPD (snStat,curPhs,nPhs,nY,nU,nP,nNodes,...
			 dvar,pvar,needF,needJ)

m = .3;
l = .5;

F = 0;
J = 0;

if (needF > 0)
  for jt = 1:nNodes,
    F(1,jt) = dvar(3,jt);
    F(2,jt) = dvar(4,jt);
    F(3,jt) = dvar(5,jt)*dvar(1,jt) / m;
    F(4,jt) = dvar(5,jt)*dvar(2,jt) / m - pvar(1);
  end
end


if (needJ > 0)
  for jt = 1:nNodes,
    J(1,jt) = 0;
    J(5,jt) = 0;
    J(9,jt) = 1;
    J(13,jt) = 0;
    J(17,jt) = 0;
    J(21,jt) = 0;
    J(25,jt) = 0;

    J(2,jt) = 0;
    J(6,jt) = 0;
    J(10,jt) = 0;
    J(14,jt) = 1;
    J(18,jt) = 0;
    J(22,jt) = 0;
    J(26,jt) = 0;

    J(3,jt) = dvar(5,jt) / m;
    J(7,jt) = 0;
    J(11,jt) = 0;
    J(15,jt) = 0;
    J(19,jt) = dvar(1,jt) / m;
    J(23,jt) = 0;
    J(27,jt) = 0;

    J(4,jt) = 0;
    J(8,jt) = dvar(5,jt) / m;
    J(12,jt) = 0;
    J(16,jt) = 0;
    J(20,jt) = dvar(2,jt) / m;
    J(24,jt) = 0;
    J(28,jt) = -1;
  end
end


