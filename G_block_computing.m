% Compute the block G^{ij} of the G matrix for ideal IAF neurons with
% delay.
%
% Inputs
%   tk1 - spike times of neuron i
%   tk2 - spike times of neuron j
%   d1  - delay of neuron i
%   d2  - delay of neuron j
%   W   - bandwidth
%   dt  - time step
%
% Output
%   G   - block G^{ij} of G matrix

% Author: Eftychios A. Pnevmatikakis
% Copyright 2009 Trustees of Columbia University

function G = G_block_computing(tk1, tk2, d1, d2, W, dt)

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
