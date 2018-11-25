% function  ivalue = sncGeti( option )
%     The  optional INTEGER-valued parameter defined by the
%     string "option" is assigned to rvalue.
%
%     For a description of the optional parameters, see
%     the snopt documentation.
%
function ivalue = sncGeti ( option )

iOpt = 7;
ivalue = ctmex ( iOpt, option );
