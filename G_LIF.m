function G = G_LIF(tk,RC)

% G_LIF create the matrix G [equation (15)]

% G = G_LIF(tk,RC) creates the matrix G with entries G[i,j] =
% <phi_i,psi_j> for the reconstruction of a stimulus that belongs in the
% L2 space and is encoded with a leaky integrate-and-fire neuron.

% Inputs
% tk:  vector of spike times
% RC:  time constant of the LIF neuron (R*C)

% Output
% G: the matrix G

% Author(s): Eftychios A. Pnevmatikakis
% Copyright 2009 Trustees of Columbia University

ln = length(tk)-1;
G = zeros(ln,ln);

etk  = exp(-diff(tk)/RC);  % exp(-(t_{k+1}-t_k)/RC)
TkTl = (ones(ln+1,1)*tk-tk'*ones(1,ln+1))/RC; % tk - tl

g=@(x)(-x.^3-6*x);
gtk = g(TkTl');

for k = 1:ln
    for l = 1:ln
        if k<l
            G(k,l) = (-gtk(l+1,k+1)+gtk(l+1,k)*etk(k)+(gtk(l,k+1)-gtk(l,k)*etk(k))*etk(l))*(RC^5);
        elseif k>l
            G(k,l) = -(-gtk(l+1,k+1)+gtk(l+1,k)*etk(k)+(gtk(l,k+1)-gtk(l,k)*etk(k))*etk(l))*(RC^5);
        else
            G(k,l) = (gtk(k+1,k)*2*etk(k)+6*(1-exp(-2*(tk(k+1)-tk(k))/RC)))*(RC^5);
        end
    end
end