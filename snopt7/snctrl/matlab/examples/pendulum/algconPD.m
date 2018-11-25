function [C,G]=algconPD (snStat,curPhs,nPhs,nC,nY,nU,nP,nNodes,...
			 dvar,pvar,needC,needG)

m = .3;
l = .5;
xf = -.231625;
yf = -.443101;

C = 0;
G = 0;

if (needC > 0)
  for jt = 1:nNodes,
    C(1,jt) = dvar(3,jt)^2 + dvar(4,jt)^2 + dvar(5,jt)*(l^2) / m - dvar(2,jt)*pvar(1);
    C(2,jt) = 0;
  end
  jt = nNodes;
  C(2,jt) = dvar(6,jt) - .5*( (dvar(1,jt)-xf)^2 + (dvar(2,jt)-yf)^2 );
end


if (needG > 0)
  for jt = 1:nNodes,
    G(1,jt) = 0;
    G(3,jt) = -pvar(1);
    G(5,jt) = 2*dvar(3,jt);
    G(7,jt) = 2*dvar(4,jt);
    G(9,jt) = l^2/m;
    G(11,jt) = 0;
    G(13,jt) = -dvar(2,jt);

    G(2,jt) = 0;
    G(4,jt) = 0;
    G(6,jt) = 0;
    G(8,jt) = 0;
    G(10,jt) = 0;
    G(12,jt) = 0;
    G(14,jt) = 0;
  end
  jt = nNodes;
  G(2,jt) = -(dvar(1,jt) - xf);
  G(4,jt) = -(dvar(2,jt) - yf);
  G(6,jt) = 0;
  G(8,jt) = 0;
  G(10,jt) = 0;
  G(12,jt) = 1;
  G(14,jt) = 0;
end


