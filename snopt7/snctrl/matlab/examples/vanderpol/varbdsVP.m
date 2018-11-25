function [lbds,ubds,x] = varbdsVP(curPhs,nPhs,nY,nU,nP,nC,nNodes)

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


% At t0
lbds(1,1) = 0;
lbds(2,1) = 0;
lbds(3,1) = 1;

ubds(1,1) = 0;
ubds(2,1) = 0;
ubds(3,1) = 1;

% At tf
lbds(2,nNodes) = 0;
lbds(3,nNodes) = -1;

ubds(2,nNodes) = 0;
ubds(3,nNodes) = -1;

% At t in (t0, tf)
for j = 2:nNodes-1,
  lbds(3,j) = -1;
  ubds(3,j) =  1;
end

for j = 1:nNodes,
  ubds(1+nY,j) = .75;
end
