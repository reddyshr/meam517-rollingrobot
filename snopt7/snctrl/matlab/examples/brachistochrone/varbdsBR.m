function [lbds,ubds,x] = varbdsBR(nPhs,curPhs,nY,nU,nP,nC,nNodes)

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

ubds(1,1) = 0;
ubds(2,1) = 0;

% At tf
lbds(1,nNodes) = bminus;
lbds(2,nNodes) = -.5;

ubds(1,nNodes) = bplus;
ubds(2,nNodes) = -.5;

