function [C,G]=algconCM (snStat,curPhs,nPhs,nC,nY,nU,nP,nNodes,...
			 dvar,pvar,needC,needG)

C = 0;
G = 0;

if (needC > 0),
  for jt = 1:nNodes,
    C(1,jt) = dvar(1,jt) + dvar(2,jt) - dvar(2+nY,jt);
  end
end


if (needG > 0),
  for jt = 1:nNodes,
    G(1,jt) =  1;
    G(2,jt) =  1;
    G(3,jt) =  0;
    G(4,jt) = -1;
  end
end


