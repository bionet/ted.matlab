function h = cross_fb(t,n1,n2,tauf,scale)

% (t,n1,n2,tauf,scale)
% calculation of the crossfeedback term (coupling function)

% cross_fb(t,n1,n2,tauf,scale) calculates the crossfeedback term for coupling 
% the two components of an ON-OFF IF neuron pair
% The mean value of the term is zero to avoid unstable behavior
% The general term is given by eq. (33)
% h = scale/tauf*exp(-t/tauf)*((t/tauf)^n1/(n1!) - (t/tauf)^n2/(n2!))

% Inputs:
% t:      time vector
% n1:     exponent for positive term
% n2:     exponent for negative term
% tauf:   time constant
% scale:  amplitude of crossfeedback term

% Output:
% h:      response of the crossfeedback term

% Author(s): Eftychios A. Pnevmatikakis
% Copyright 2010 Trustees of Columbia University

h = scale/tauf*exp(-t/tauf).*((t/tauf).^n1/factorial(n1)-(t/tauf).^n2/factorial(n2));