%G_POP_IF Compute the reconstruction matrix G for consistent recovery.
%   G = G_POP_IF(TK,LN,W) computes the matrix G with entries 
%   G[i,j] = <phi_k^i,psi_l^j> used to reconstruct a signal in L2
%   space that was encoded with a population of ideal IAF
%   neurons. TK contains the spike times, while LN denotes the
%   number of spikes due to each neuron. If specified, W denotes
%   the weights on the inputs prior to encoding; if not specified,
%   the weights are all assumed to be 1.
%
%   The calculation is described in further detail in Equation 29 
%   of the Consistent Recovery paper mentioned in the toolbox 
%   references.

%   Author: Eftychios A. Pnevmatikakis
%   Copyright 2009-2010 Trustees of Columbia University

function G = G_pop_IF(TK,ln,varargin)

N = length(ln);
G = zeros(sum(ln),sum(ln));
ln2 = cumsum([0,ln]);

if nargin > 2
    w = varargin{1};
else
    w = ones(1,N);
end

for i = 1:N
    for j = 1:N
        Gb = G_block_IF(TK(1:ln(i)+1,i)',TK(1:ln(j)+1,j)');
        G(ln2(i)+1:ln2(i+1),ln2(j)+1:ln2(j+1)) = w(i)*w(j)*Gb;
    end
end
