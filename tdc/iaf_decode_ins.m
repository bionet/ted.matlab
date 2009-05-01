%IAF_DECODE_INS Decode a signal encoded with an IAF neuron
%   U_REC = IAF_DECODE_INS(S,T,W,B,R,C) decodes the signal encoded as
%   a binary vector S defined over the times T using an IAF neuron
%   with bias B, resistance R, and capacitance C. The intervals
%   between the successive elements of T are assumed to be
%   uniform. The original signal's bandwidth W (in rad/s) must also be
%   specified. Note that the decoding algorithm does not need the
%   threshold of the encoder to recover the encoded signal.
%
%   U_REC = IAF_DECODE_INS(S,T,W,B,C) decodes the encoded signal S
%   assuming a non-leaky IAF neuron model (i.e., R = inf). If C is
%   not specified, it is assumed to be 1.

%   Author(s):     L. Givon, E. Pnevmatikakis
%   Copyright 2008 L. Givon, E. Pnevmatikakis
%                  Bionet Research Group
%                  Dept. of Electrical Engineering
%                  Columbia University

function u_rec = iaf_decode_ins(s,t,W,b,varargin)

R = inf;
if nargin == 4,
  C = 1;
elseif nargin == 5,
  C = varargin{1};
elseif nargin == 6,
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
G = zeros(nsk-1,nsk-1);
for i=1:nsk-1,
  for j=1:nsk-1,
    G(i,j) = dt*sum(sinc(Wpi*(t(tk(i):tk(i+1)) - t(sk(j))))... 
                    *Wpi.*exp(-(t(tk(i+1))-t(tk(i):tk(i+1)))./RC)); 
  end
end
    
% Apply compensation principle:
B = diag(ones(1,nsk-1),1)-eye(nsk);
B_inv = inv(B);
if isinf(R),
  Bq = -b*diff(diff(t(tk)))';
else
  Bq = RC*b*diff(exp(-diff(t(tk))/RC))';
end

% Discard the last rows and columns of B and B_inv to
% eliminate boundary issues:
B = B(1:end-1,1:end-1);
B_inv = B_inv(1:end-1,1:end-1);

% Reconstruct signal by adding up the weighted sinc functions:
u_rec = zeros(1,nt);
c = B_inv*pinv(B*G*B_inv)*Bq;
for i=1:nsk-1,
  u_rec = u_rec + sinc(Wpi*(t-t(sk(i))))*Wpi*c(i);
end
