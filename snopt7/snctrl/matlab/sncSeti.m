% function sncSeti( option, ivalue )
%     Sets the INTEGER-VALUED optional parameter defined by the
%     string "option" is assigned to rvalue.
%
%     For a description of all the optional parameters, see the
%     snopt documentation.
%
%     Do not try to set the unit number of the summary or print file.
%     Use the MATLAB functions snsummary and snprintfile instead.
%
function sncSeti( option, ivalue )

iOpt = 3;
snctrlmex( iOpt, option, ivalue );
