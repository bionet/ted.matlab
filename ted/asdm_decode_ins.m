%ASDM_DECODE_INS Decode a signal encoded with an asyncrhonous sigma-delta modulator.
%   U_REC = ASDM_DECODE(S,DUR,DT,BW,B) decodes the signal with
%   duration DUR s and bandwidth BW rad/s encoded as an array of
%   spike intervals S using an asynchronous sigma-delta modulator 
%   with bias B. The recovered signal is assumed to be sampled
%   at sampling rate 1/DT Hz.
%
%   Note that the decoding algorithm does not need the 
%   threshold and integration constant of the encoder to recover
%   the encoded signal.

%   Author(s): Lev Givon
%   Copyright 2009 Trustees of Columbia University

function u_rec = asdm_decode_ins(s,dur,dt,bw,b)

ns = length(s);
if ns < 2,
  error('s must contain at least 2 elements');
end

ts = cumsum(s);
tsh = (ts(1:end-1)+ts(2:end))/2;

nt = floor(dur/dt);
t = linspace(0,dur,nt);
nsh = length(tsh);

bwpi = bw/pi;

% Compute G matrix:
G = zeros(nsh,nsh);
for j=1:nsh,
  temp = si(bw*(ts-tsh(j)))/pi;
  for i=1:nsh,
    G(i,j) = temp(i+1)-temp(i);
  end
end

% Apply compensation principle:
B = diag(ones(1,nsh-1),-1)+eye(nsh);
Bq = ((-1).^[0:nsh-1].*(s(2:end)-s(1:end-1))*b)';

% Reconstruct signal by adding up the weighted sinc functions; the
% first row of B is removed to eliminate boundary issues. The
% weighted sinc functions are computed on the fly here to save memory:
u_rec = zeros(1,nt);
c = pinv(B(2:end,:)*G)*Bq(2:end);
for i=1:nsh,
  u_rec = u_rec + sinc(bwpi*(t-tsh(i)))*bwpi*c(i);
end
