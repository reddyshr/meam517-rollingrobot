function runVanderpol
%
%
%

name = 'vanderpl';
nY = 3;
nU = 1;
nP = 0;
nC = 0;
objL(1) = 1;

nPhs = 1;
npInt(1) = 10;
phsPt(1) = 0;
phsPt(2) = 5;

ytype(1,1) = 0;
ytype(2,1) = 0;
ytype(3,1) = 1;

ctype = 0;


% Initialize and set print file.
sncInit ( 'vanderpol.out' );

% Read a specs file.
specFile.spc = which('vanderpol.spc');
sncSpec ( specFile.spc );


% Solve the problem.
Start = 0;
snctrlD ( Start, name, nY, nU, nP, nC, nPhs, objL, npInt, phsPt, ...
	  ytype, ctype, 'varbdsVP', 'odeconVP', '' );
