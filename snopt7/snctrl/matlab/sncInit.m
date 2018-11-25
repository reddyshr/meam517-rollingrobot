function sncInit ( file )
%
% This is the initialization routine for the
% optimal control interface snctrl.  The input string
% 'file' should define the name of the print file.
%

iOpt = 1;
ctmex ( iOpt, file );
