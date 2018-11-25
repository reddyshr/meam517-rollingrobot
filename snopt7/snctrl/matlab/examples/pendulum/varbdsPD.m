function [lbds,ubds,x,plbds,pubds,p,clbds,cubds] = varbdsPD(curPhs,nPhs,nY,nU,nP,nC,nNodes)

bminus = -10^20;
bplus  =  10^20;

xf = -.231625;
yf = -.443101;
l  = .5;

% Default set everything to free bounds and initial point is 0.
for j = 1:nNodes,
  for i = 1:nY+nU,
    lbds(i,j) = bminus;
    ubds(i,j) = bplus;
    x(i,j)    = 0;
  end
end

for j = 1:nNodes,
  x(1,j) = .4;
  x(2,j) = -.3;
  x(3,j) = 0;
  x(4,j) =0;
  x(5,j) = -5;
  x(6,j) =0;
end
x(6,nNodes) = .5*( (x(1,nNodes)-xf)^2 + (x(2,nNodes)-xf)^2 );


% At t0
lbds(1,1) = .4;
lbds(2,1) = -.3;
lbds(3,1) = 0;
lbds(4,1) = 0;

ubds(1,1) = .4;
ubds(2,1) = -.3;
ubds(3,1) = 0;
ubds(4,1) = 0;

% Parameters
plbds(1) = .01;
pubds(1) = 100;
p(1) = 20;

% Algebraic constraint
clbds(1,1) = 0;
cubds(1,1) = 0;

clbds(2,1) = 0;
cubds(2,1) = 0;
