%ASDM_RECOVERABLE Check whether TDM parameters permit signal recovery.
%   ASDM_RECOVERABLE(U,BW,B,D,K) returns TRUE if the signal U can
%   be recovered using an ASDM decoder with capacitance K,
%   threshold D, bias B, and bandwidth BW (in rad/s).
%   Note that the bound assumed by this check is not as strict as
%   that described in most of Prof. Lazar's papers.

%   Author: Lev Givon
%   Copyright 2009-2010 Trustees of Columbia University

function result = asdm_recoverable(u,bw,b,d,k)

result = false;
c = max(abs(u));
if c >= b,
  fprintf(1,'bias too low\n');
elseif 2*k*d/(b)*bw/pi >= 1,
  fprintf(1,'reconstruction condition not satisfied\n');
  fprintf(1,'try increasing b or reducing k or d\n');
else
  result = true;
end
