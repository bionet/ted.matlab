%IAF_DECODE_FAST Decode a signal encoded with an IAF neuron.
%   U_REC = IAF_DECODE_FAST(S,DUR,DT,BW,M,B,D,R,C) decodes the signal
%   with duration DUR s and bandwidth BW rad/s encoded as an array
%   of spike intervals S using an integrate-and-fire neuron with
%   bias B, firing threshold D, resistance R, and capacitance
%   C. The recovered signal is assumed to be sampled sampling rate
%   1/DT Hz. This function uses a fast decoding algorithm; M
%   specifies the number of bins used by this algorithm.
%
%   The parameters R and C are optional. If R = inf (the default), an
%   ideal neuron model is used. C is assumed to be 1 if not
%   specified.

%   Author: Lev Givon
%   Copyright 2009-2015 Lev Givon

function u_rec = iaf_decode_fast(s,dur,dt,bw,M,b,d,varargin)

R = inf;
C = 1;
if nargin >= 8,
  R = varargin{1};
end
if nargin >= 9,
  R = varargin{1};
  C = varargin{2};
end
if nargin >= 10,
  error('Too many input arguments.');
end

ns = length(s);
if ns < 2,
  error('s must contain at least 2 elements.');
end

ts = cumsum(s);
tsh = (ts(1:end-1)+ts(2:end))/2;
nsh = length(tsh);

t = [0:dt:dur];

RC = R*C;
jbwM = j*bw/M;

% Compute quanta:
if isinf(R),
  q = (C*d - b*s(2:end))';
else
  q = C*(d + b*R*(exp(-s(2:end)/RC)-1))';
end

% Compute approximation coefficients:
a = bw/(pi*(2*M+1));
m = [-M:M];
P_inv = -triu(ones(nsh));
S = exp(-jbwM*m'*ts(1:end-1));
D = diag(s(2:end));
SD = S*D;
T = a*SD*S';
dd = a*pinv(T)*SD*P_inv*q;

% Reconstruct signal: 
u_rec = real(jbwM*(m.*dd.')*exp(jbwM*m'*t));
