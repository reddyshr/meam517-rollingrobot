function [lbds,ubds,x,plbds,pubds,p,clbds,cubds] = varbdsVT(curPhs,nPhs,nY,nU,nP,nC,nNodes)

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
lbds(3,1) = 0;
lbds(4,1) = 0;

ubds(1,1) = 0;
ubds(2,1) = 0;
ubds(3,1) = 0;
ubds(4,1) = 0;

% At tf
lbds(1,nNodes) = 1;
lbds(2,nNodes) = 0;

ubds(1,nNodes) = 1;
ubds(2,nNodes) = 2;

% Parameters
plbds(1) = 0;
pubds(1) = bplus;
p(1)     = 1;

% Algebraic constraint
clbds(1,1) = bminus;
cubds(1,1) = .2;

