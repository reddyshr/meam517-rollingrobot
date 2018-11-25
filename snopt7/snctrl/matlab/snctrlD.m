function [INFO,x,hs]=snctrlD( Start,name,nY,nU,nP,ytype,odecon,varbds, ...
			      nPhs,npInt,phsPt,objL,varargin )
%function [INFO,x,hs]=snctrlD( Start,name,nY,nU,nP,ytype,odecon,varbds, ...
%			      nPhs,npInt,phsPt,objL,varargin)
%
% [INFO,x,hs]=snctrlD( Start,name,nY,nU,nP,ytype,odecon,varbds, ...
%     	               nPhs,npInt,phsPt,objL)
% [INFO,x,hs]=snctrlD( Start,name,nY,nU,nP,ytype,odecon,varbds,...
%		       nPhs,npInt,phsPt,objL,nC,ctype,algcon)
%
%

iOpt = 1;
if nargin == 12,
  [INFO,x,hs]=snctrlmex( iOpt,Start,name,nY,nU,nP,ytype,odecon,varbds, ...
			 nPhs,npInt,phsPt,objL )

elseif nargin == 15
  [INFO,x,hs]=snctrlmex( iOpt,Start,name,nY,nU,nP,ytype,odecon,varbds, ...
			 nPhs,npInt,phsPt,objL,nC,ctype,algcon )

end