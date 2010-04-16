%IAF_DECODE Decode a signal encoded with an IAF neuron.
%   U_REC = IAF_DECODE(S,DUR,DT,BW,B,D,R,C) decodes the signal with
%   duration DUR s and bandwidth BW rad/s encoded as an array of
%   spike intervals S using an integrate-and-fire neuron with bias
%   B, firing threshold D, resistance R, and capacitance C. The
%   recovered signal is assumed to be sampled at sampling rate 1/DT
%   Hz. 
%
%   The parameters R and C are optional. If R = inf (the default), an
%   ideal neuron model is used. C is assumed to be 1 if not
%   specified.

%   Author: Lev Givon
%   Copyright 2009-2010 Trustees of Columbia University

function u_rec = iaf_decode(s,dur,dt,bw,b,d,varargin)

R = inf;
C = 1;
if nargin >= 7,
  R = varargin{1};
end
if nargin >= 8,
  R = varargin{1};
  C = varargin{2};
end
if nargin >= 9,
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

bwpi = bw/pi;
RC = R*C;

% Compute G matrix and quanta:
G = zeros(nsh,nsh);
if isinf(R),
  for j=1:nsh,
    temp = si(bw*(ts-tsh(j)))/pi;
    for i=1:nsh,
      G(i,j) = temp(i+1)-temp(i);
    end
  end
  q = (C*d-b*s(2:end))';
else
  for i=1:nsh,
    for j=1:nsh,
      G(i,j) = quad(@(t) sinc(bwpi*(t-tsh(j)))*bwpi.*exp((ts(i+1)-t)/-RC),...
                    ts(i),ts(i+1));
    end
  end
  q = (C*(d+b*R*(exp(-s(2:end)/RC)-1)))';
end
G_inv = pinv(G);

% Reconstruct signal by adding up the weighted sinc functions:
u_rec = zeros(1,length(t));
c = G_inv*q;
for i=1:nsh,
  u_rec = u_rec + sinc(bwpi*(t-tsh(i)))*bwpi*c(i);
end

