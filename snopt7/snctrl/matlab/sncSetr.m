% function sncSetr( option, rvalue )
%     Sets the optional REAL-VALUED parameter defined by the
%     string "option" to the value  rvalue.
%
%     For a description of all the optional parameters, see the
%     snopt documentation.
%
%     Do not try to set the unit number of the summary or print file.
%     Use the MATLAB functions snsummary and snprintfile instead.
%
function sncSetr( option, rvalue )

iOpt = 4;
snctrlmex( iOpt, option, rvalue );
