function p = p_IF(TK,RC,varargin)

% p_IF creates the vector of inner product measurements p_k = <1,\phi_k>

% p_IF(TK,RC) creates the vector for a LIF neuron with time constant RC
% If the neuron is an ideal integrator then choose RC = Inf
% p_IF(TK,RC,ln) creates the vector for a population of length(ln) neurons 
% where ln represents the number of spikes from each neuron
% p_IF(TK,RC,ln, w) creates the vector for a population of length(ln) neurons 
% where ln represents the number of measurements from each neuron, and each
% neuron j weights the input by a number w(j) before encoding

% Inputs
% TK:   vector of spike trains
% RC:   time constants of each neuron
% ln:   # of measurements for each neuron
% w:    weighing factor for each neuron

% Output
% p:    the p-vector

% Author(s): Eftychios A. Pnevmatikakis
% Copyright 2010 Trustees of Columbia University

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