function G = G_pop_IF(TK,ln,varargin)

% G_pop_IF create the matrix G [equation (29)]

% G = G_pop_IF(tk,ln) creates the matrix G with entries G[i,j] =
% <phi_k^i,psi_l^j> for the reconstruction of a stimulus that belongs in the
% L2 space and is encoded with a population of ideal IF neurons

% G = G_pop_IF(tk,ln,w) accounts for any weighting of the inputs prior to
% encoding

% Inputs
% TK:  vector of spike times from all neurons
% ln:  number of measurements from each neuron (#spikes-1)
% w:   weighing vector

% Output
% G: the matrix G

% Author(s): Eftychios A. Pnevmatikakis
% Copyright 2009 Trustees of Columbia University

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
        Gb = Gblock_IF(TK(1:ln(i)+1,i)',TK(1:ln(j)+1,j)');
        G(ln2(i)+1:ln2(i+1),ln2(j)+1:ln2(j+1)) = w(i)*w(j)*Gb;
    end
end