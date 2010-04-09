%G_IF Create the reconstruction matrix G.
%   G = G_IF(TK,RC) computes the reconstruction matrix G with entries
%   G[i,j] = <phi_i,psi_j> used to decode a signal in L2 space that
%   was encoded by a leaky IAF neuron as the spike times TK with time
%   constant RC; if no time constant is specified, the neuron is
%   assumed to be ideal.
%   
%   The calculation is described in further detail in Equations
%   17 and 19 of the Consistent Recover paper mentioned in the
%   toolbox references.

%   Author: Eftychios A. Pnevmatikakis
%   Copyright 2009-2010 Trustees of Columbia University

function G = G_IF(tk,varargin)

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