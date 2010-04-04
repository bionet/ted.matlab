%P_IF Compute the inner product <1,\phi_k>.
%   P = P_IF(TK,RC) computes the inner product <1,\phi_k> for
%   a leaky IAF neuron with time constant RC (an ideal neuron may
%   be simulated by setting RC to Inf). The spike times generated
%   by the neuron are specified in TK.
%
%   P = P_IF(TK,RC,LN,W) computes the inner product for a population
%   of IAF neurons; LN contains the number of spikes from each
%   neuron, and W contains the values by which each neuron weights
%   the input. If W is not specified, the input is not weighted.

%   Author: Eftychios A. Pnevmatikakis
%   Copyright 2009-2010 Trustees of Columbia University

function p = p_IF(TK,RC,varargin)

if nargin > 2
    ln = varargin{1};
else
    ln = length(TK)-1;
    w = 1;
    N = 1;
end
if nargin > 3
    w = varargin{2};
    N = length(ln);
else
    w = ones(1,length(ln));
    N = length(ln);
end
if nargin > 4
    error('Too many input arguments.');
end

ln2 = cumsum([0,ln]);
p = zeros(sum(ln),1);

if isinf(RC)
    for i = 1:N
        p(ln2(i)+1:ln2(i+1)) = w(i)*diff(TK(1:ln(i)+1,i));
    end
else
    for i = 1:N
        p(ln2(i)+1:ln2(i+1)) = w(i)*RC(i)*(1-exp(-diff(TK(1:ln(i)+1,i))/RC(i)));
    end
end