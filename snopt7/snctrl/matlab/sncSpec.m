% function inform = sncSpec( filename )
%     Causes snopt to read its optional parameters from the named file.
%     The format of this file is described in the snopt documentation.
%
%     Returns 0 if successful, and a positive number otherwise.
function INFO = sncSpec( buffer )

iOpt = 9;
INFO = snctrlmex( iOpt, buffer );
