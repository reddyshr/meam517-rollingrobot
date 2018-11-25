function runCatmix
%
%
%

name = 'catmix  ';
nY = 2;
nU = 2;
nP = 0;
nC = 1;
objL(1) = -2;

nPhs = 1;
npInt(1) = 10;
phsPt(1) = 0;
phsPt(2) = 1;

ytype(1,1) = 0;
ytype(2,1) = 0;

ctype(1,1) = 1;

% Initialize the control interface.  No print file.
sncInit ('');

% Solve the problem.
Start = 0;
snctrlD ( Start, name, nY, nU, nP, nC, nPhs, objL, ...
	  npInt, phsPt, ytype, ctype, 'varbdsCM', 'odeconCM', 'algconCM' );
