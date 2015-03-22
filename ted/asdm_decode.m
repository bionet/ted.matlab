%ASDM_DECODE Decode a signal encoded with an asynchronous sigma-delta modulator.
%   U_REC = ASDM_DECODE(S,DUR,DT,BW,B,D,K) decodes the signal with
%   duration DUR s and bandwidth BW rad/s encoded as an array of spike
%   intervals S using an asynchronous sigma-delta modulator with bias
%   B, Schmitt trigger threshold D, and integration constant K. The
%   recovered signal is assumed to be sampled at sampling rate 1/DT Hz.

%   Author: Lev Givon
%   Copyright 2009-2015 Lev Givon

function u_rec = asdm_decode(s,dur,dt,bw,b,d,k)

ns = length(s);
if ns < 2,
  error('s must contain at least 2 elements');
end

ts = cumsum(s);
tsh = (ts(1:end-1)+ts(2:end))/2;
nsh = length(tsh);

t = [0:dt:dur];

bwpi = bw/pi;

% Compute G matrix:
G = zeros(nsh,nsh);
for j=1:nsh,
  temp = si(bw*(ts-tsh(j)))/pi;
  G(:,j) = temp(2:end)-temp(1:end-1);
end
G_inv = pinv(G);

% Compute quanta:
q = ((-1).^[1:nsh].*(2*k*d-b*s(2:end)))';

% Reconstruct signal by adding up the weighted sinc functions. The
% weighted sinc functions are computed on the fly here to save memory:
u_rec = zeros(1,length(t));
c = G_inv*q;
for i=1:nsh,
  u_rec = u_rec + sinc(bwpi*(t-tsh(i)))*bwpi*c(i);
end
