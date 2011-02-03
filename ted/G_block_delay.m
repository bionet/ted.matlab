%G_BLOCK_DELAY Compute part of the reconstruction matrix for delayed IAF neurons.
%   G = COMPUTE_G_BLOCK(TK1,TK2,D1,D2,W,DT) computes the block G^{ij}
%   of the reconstruction matrix used to decode signals with bandwidth
%   W rad/s sampled at a rate 1/DT Hz encoded by ideal IAF neurons
%   with delays. TK1 and TK2 contain the spike times of neurons i and
%   j, respectively; D1 and D2 contain the delays of neurons i and j.

%   Author: Eftychios A. Pnevmatikakis
%   Copyright 2009-2011 Eftychios A. Pnevmatikakis

function G = G_block_delay(tk1,tk2,d1,d2,W,dt)

if nargin ~= 6,
  error('Wrong number of input arguments.');
end

sk2 = (tk2(1:end-1) + tk2(2:end))./2; % midpoints of interspike intervals

G = zeros(length(tk1)-1, length(sk2));
for k=1:length(tk1)-1
    for l=1:length(sk2)
        G(k,l) = dt*trapz(sinc(W/pi*((tk1(k):dt:tk1(k+1))-sk2(l)-d1+d2)))*W/pi;
    end
end
