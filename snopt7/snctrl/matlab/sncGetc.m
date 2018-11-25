% function  cvalue = sncGetC( option )
%     The  optional parameter defined by the
%     string "option" is assigned to cvalue.
%
%     For a description of the optional parameters, see
%     the snopt documentation.
%
function cvalue = sncGetc ( option )

iOpt = 6;
cvalue = snctrlmex( iOpt, option );
