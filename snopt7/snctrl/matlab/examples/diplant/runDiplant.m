function runDiplant
%
%
%

name = 'diplant ';
nY = 3;
nU = 1;
nP = 1;
nC = 0;
objL(1) = 3;

nPhs = 1;
npInt(1) = 10;
phsPt(1) = 0;
phsPt(2) = 1;

ytype(1,1) = 0;
ytype(2,1) = 0;
ytype(3,1) = 1;

ctype = 0;

% Initialize and set the print file.
sncInit ( 'diplant.out' );

% Read in a specs file.
specFile.spc = which('diplant.spc');
sncSpec ( specFile.spc );

% Solve the problem.
Start = 0;
snctrlD ( Start, name, nY, nU, nP, nC, nPhs, objL, npInt, phsPt, ...
	  ytype, ctype, 'varbdsDI', 'odeconDI', '' );
