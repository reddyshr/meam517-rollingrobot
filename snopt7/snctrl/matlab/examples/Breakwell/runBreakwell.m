function runBreakwell
%
%
%

name = 'Breakwel';
nY = 3;
nU = 1;
nP = 0;
nC = 0;
objL(1) = 1;

nPhs = 1;
npInt(1) = 10;
phsPt(1) = 0;
phsPt(2) = 1;

ytype(1,1) = 0;
ytype(2,1) = 1;
ytype(3,1) = 1;

ctype = 0;

% First call sncInit to initialize and set the print file.
sncInit ('Breakwell.out');

% Read a Specs file.
specFile.spc = which('Breakwell.spc');
sncSpec ( specFile.spc );

% Now solve the proble by calling snctrl.
Start = 0;
[x, hs] = snctrlD ( Start, name, nY, nU, nP, nC, nPhs, objL, ...
		    npInt, phsPt, ytype, ctype, ...
		    'varbdsBW', 'odeconBW', '' );

