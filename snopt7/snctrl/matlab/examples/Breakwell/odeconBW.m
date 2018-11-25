function [F,J]=odeconBW (snStat,curPhs,nPhs,nY,nU,nP,nNodes,...
			 dvar,pvar,needF,needJ)

F = 0;
J = 0;

if (needF > 0)
  for jt = 1: nNodes,
    F(1,jt) = ( dvar(1+nY,jt)^2 ) / 2;
    F(2,jt) = dvar(3,jt);
    F(3,jt) = dvar(1+nY,jt);
  end
end


if (needJ > 0)
  for jt = 1:nNodes,
    J(1,jt) = 0;
    J(4,jt) = 0;
    J(7,jt) = 0;
    J(10,jt) = dvar(1+nY,jt);

    J(2,jt) = 0;
    J(5,jt) = 0;
    J(8,jt) = 1;
    J(11,jt) = 0;

    J(3,jt) = 0;
    J(6,jt) = 0;
    J(9,jt) = 0;
    J(12,jt) = 1;
  end

end


