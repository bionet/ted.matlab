%IAF_RECOVERABLE Check whether TDM parameters permit signal recovery
%   IAF_RECOVERABLE(U,W,B,D,R,C) returns TRUE if the signal U can be
%   recovered using an IAF decoder with bias B, threshold D,
%   resistance R, capacitance C, and bandwidth W (in rad/s).
%
%   Author: Lev Givon
%   Copyright 2009-2012 Lev Givon

function res = iaf_recoverable(u,W,b,d,R,C)

res = 0;
c = max(abs(u));
if c >= b,
  fprintf(1,'bias too low\n');
end
r = R*C*log(1-d/(d-(b-c)*R))*W/pi;
e = d/((b-c)*R);
if ~isreal(r),
  fprintf(1,'reconstruction condition not satisfied\n');
  fprintf(1,'r = %f + %fj\n',real(r),imag(r));
elseif r >= (1-e)/(1+e),
  fprintf(1,'reconstruction condition not satisfied\n');
  fprintf(1,'r = %f, (1-e)/(1+e) = %f\n',r,(1-e)/(1+e));
else
  res = 1;
end
