%ASDM_DECODE_FAST Decode a signal encoded with an asynchronous sigma-delta modulator.
%   U_REC = ASDM_DECODE_FAST(S,DUR,DT,BW,M,B,D,K) decodes the signal
%   with duration DUR s and bandwidth BW rad/s encoded as an array
%   of spike intervals S using an asynchronous sigma-delta modulator with
%   bias B, firing threshold D, and integration constant
%   K. The recovered signal is assumed to be sampled sampling rate
%   1/DT Hz. This function uses a fast decoding algorithm; M
%   specifies the number of bins used by this algorithm.

%   Author: Lev Givon
%   Copyright 2009-2011 Lev Givon

function u_rec = asdm_decode_fast(s,dur,dt,bw,M,b,d,k)

ns = length(s);
if ns < 2,
  error('s must contain at least 2 elements.');
end

ts = cumsum(s);
tsh = (ts(1:end-1)+ts(2:end))/2;
nsh = length(tsh);

t = [0:dt:dur];

jbwM = j*bw/M;

% Compute quanta:
q = ((-1).^[1:nsh].*(2*k*d-b*s(2:end)))';

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
