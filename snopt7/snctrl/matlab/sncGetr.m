% function  rvalue = sncGetr( option )
%     The REAL-VALUED optional parameter defined by the string "option"
%     is assigned to rvalue.
%
%     For a description of all the optional parameters, see the
%     snopt documentation.
%
function rvalue = sncGetr ( option )

iOpt = 8;
rvalue = snctrlmex( iOpt, option );
