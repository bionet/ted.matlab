function G = G_IF(tk,varargin)

% G_IF create the matrix G [equation (15)]

% G = G_IF(tk) creates the matrix G with entries G[i,j] =
% <phi_i,psi_j> for the reconstruction of a stimulus that belongs in the
% L2 space and is encoded with an ideal integrate-and-fire neuron (equation (19)).

% G = G_IF(tk,RC) creates the matrix G with entries G[i,j] =
% <phi_i,psi_j> for the reconstruction of a stimulus that belongs in the
% L2 space and is encoded with a leaky integrate-and-fire neuron with time 
% constant RC (equation (17)).

% Inputs
% tk:  vector of spike times
% RC:  time constant of the LIF neuron (R*C)

% Output
% G: the matrix G

% Author(s): Eftychios A. Pnevmatikakis
% Copyright 2009 Trustees of Columbia University

ln = length(tk)-1;
G = zeros(ln,ln);

if nargin > 1
    RC = varargin{1};
else
    RC = Inf;
end

if isinf(RC)
    for i=1:ln
        for j=1:ln
            tmz=tk(min(i,j));
            tmp=tk(min(i,j)+1);
            tpz=tk(max(i,j));
            tpp=tk(max(i,j)+1);
            G(i,j)=((tpp-tmz)^5-(tpp-tmp)^5+(tpz-tmp)^5-(tpz-tmz)^5)/20;
        end
    end    
else
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
end