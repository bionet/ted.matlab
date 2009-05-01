%IAF_DECODE_REC Decode a signal encoded with an IAF neuron
%   U_REC = IAF_DECODE(S,T,W,L,B,D,R,C) decodes the signal encoded as
%   a binary vector S defined over the times T using an IAF
%   neuron. The decoding is performed assuming the IAF bias B,
%   threshold D, resistance R, and capacitance C. The intervals
%   between the successive elements of T are assumed to be
%   uniform. The original signal's bandwidth W (in rad/s) must also be
%   specified. L indicates the number of times the time decoding
%   operator must be recursively applied to reconstruct the signal.
%
%   U_REC = IAF_DECODE(S,T,W,L,B,D,C) decodes the encoded signal S
%   assuming a non-leaky IAF neuron model (i.e., R = inf). If C is
%   not specified, it is assumed to be 1.

%   Author(s):     L. Givon, E. Pnevmatikakis, R. Turetsky
%   Copyright 2008 L. Givon, E. Pnevmatikakis, R. Turetsky
%                  Bionet Research Group
%                  Dept. of Electrical Engineering
%                  Columbia University

function u_rec = iaf_decode_rec(s,t,W,L,b,d,varargin)

R = inf;
if nargin == 6,
  C = 1;
elseif nargin == 7,
  C = varargin{1};
elseif nargin == 8,
  R = varargin{1};
  C = varargin{2};
else
  error('too many arguments');
end

nt = length(t);
if nt < 2,
  error('t must contain at least 2 elements');
end

tk = find(s);
sk = round((tk(1:end-1) + tk(2:end))./2);
dt = t(2)-t(1);

nsk = length(sk);

Wpi = W/pi;
RC = R*C;

% Compute shifted sincs and pseudoinverse matrix:
g = zeros(nt,nsk);
G = zeros(nsk,nsk);
for i=1:nsk,
  g(:,i) = sinc(Wpi*(t-t(sk(i))))'*Wpi;
  for j=1:nsk,
    G(i,j) = dt*sum(sinc(Wpi*(t(tk(i):tk(i+1)) - t(sk(j))))...
                    *Wpi.*exp(-(t(tk(i+1))-t(tk(i):tk(i+1)))./RC));
  end
end
IG = eye(nsk)-G;

% Compute quanta:
if isinf(R),
  q = C*d - b*diff(t(tk))';
else
  q = C*(d + b*R*(exp(-diff(t(tk))'/RC)-1));
end

% Recursively reconstruct signal: 
u_rec = zeros(1,nt);
IGj = eye(nsk);
for j=0:L,
  u_curr = (g*IGj*q)';
  if j>0,
    u_rec = u_curr+u_rec;
  else,
    u_rec = u_curr;
  end
  IGj = IGj*IG;
end

