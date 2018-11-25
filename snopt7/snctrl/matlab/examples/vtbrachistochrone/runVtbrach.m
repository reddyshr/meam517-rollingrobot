function runVtbrach
%
%
%

name = 'vtbrach ';
nY = 4;
nU = 1;
nP = 1;
nC = 1;
objL(1) = 4;

nPhs = 1;
npInt(1) = 10;
phsPt(1) = 0;
phsPt(2) = 1;

ytype(1,1) = 0;
ytype(2,1) = 0;
ytype(3,1) = 0;
ytype(4,1) = 1;

ctype(1,1) = 1;

% Initialize and set the print file.
sncInit ( 'vtbrach.out' );

% Read specs file.
specFile.spc = which('vtbrach.spc');
sncSpec ( specFile.spc );

% Solve the problem.
Start = 0;
snctrlD ( Start, name, nY, nU, nP, nC, nPhs, objL, npInt, phsPt, ...
	  ytype, ctype, 'varbdsVT', 'odeconVT', 'algconVT' );
