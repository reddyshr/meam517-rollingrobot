function runBrachistochrone
%
%
%

name = 'brachist';
nY = 2;
nU = 1;
nP = 0;
nC = 0;

objL(1) = 1;
%objL(2) = 1;

nPhs = 1;
npInt(1) = 10;
%npInt(2) = 5;

phsPt(1) = 0;
%phsPt(2) = .5;
phsPt(2) = 1;


ytype(1,1) = 0;
ytype(2,1) = 1;

%ytype(1,2) = 0;
%ytype(2,2) = 1;

ctype = 0;

% First call sncInit and set the print file.
sncInit ('brach.out');

% Set to HS discretization and Print Control Summary
sncSet ( 'Discretization HS');
sncSet ( 'Control Solution Yes');


% Now solve the proble by calling snctrl.
Start = 0;
snctrlD ( Start, name, nY, nU, nP, nC, nPhs, objL, npInt, phsPt, ...
	  ytype, ctype, 'varbdsBR', 'odeconBR', '' );

