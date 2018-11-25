function runPendulum
%
%
%

name = 'Pendulum';
nY = 4;
nU = 2;
nP = 1;
nC = 2;
objL(1) = -2;

nPhs = 1;
npInt(1) = 25;
phsPt(1) = 0;
phsPt(2) = 2;

ytype(1,1) = 1;
ytype(2,1) = 1;
ytype(3,1) = 0;
ytype(4,1) = 0;

ctype(1,1) = 0;
ctype(2,1) = 0;

% Initialize the interface and set the print file.
sncInit ( 'pendulum.out' );

% Read a specs file.
specFile.spc = which('pendulum.spc');
sncSpec ( specFile.spc );

% Solve the problem.
Start = 0;
snctrlD ( Start, name, nY, nU, nP, nC, nPhs, objL, npInt, phsPt, ...
	  ytype, ctype, 'varbdsPD', 'odeconPD', 'algconPD' );
