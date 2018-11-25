% function  ivalue = sncGet( option )
%     The  optional INTEGER-valued parameter defined by the
%     string "option" is assigned to ivalue.
%
%     For a description of the optional parameters, see
%     the snopt documentation.
%
function ivalue = sncGet( option )

iOpt = 5;
ivalue = snctrlmex( iOpt, option );
