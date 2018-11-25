function [lbds,ubds,x,clbds,cubds] = varbdsCM(curPhs,nPhs,nY,nU,nP,nC,nNodes)

bminus = -10^20;
bplus  =  10^20;

% Default set everything to free bounds and initial point is 0.
for j = 1:nNodes,
  for i = 1:nY+nU,
    lbds(i,j) = bminus;
    ubds(i,j) = bplus;
    x(i,j)    = 0;
  end
end


% At t in (t0, tf)
for j = 1:nNodes,
  lbds(1+nY,j) = 0;
  ubds(1+nY,j) = 1;
end

% At t0
lbds(1,1) = 1;
lbds(2,1) = 0;

ubds(1,1) = 1;
ubds(2,1) = 0;

% Algebraic constraint
clbds(1,1) = 1;
cubds(1,1) = 1;