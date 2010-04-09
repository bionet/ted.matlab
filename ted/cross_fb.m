%CROSS_FB Compute cross-feedback for ON-OFF IAF neuron pair.
%   H = CROSS_FB(T,N1,N2,TAUF,SCALE) computes the cross-feedback
%   term coupling the two components of an ON-OFF IAF neuron
%   pair over the times T. N1 and N2 denote the exponents for the
%   positive and negative feedback terms, respectively; TAUF
%   denotes the time constant, and SCALE the amplitude of the
%   cross-feedback.
%
%   The calculation is described in further detail in Equation 33 of
%   the Consistent Recovery paper mentioned in the toolbox references.

%   Author: Eftychios A. Pnevmatikakis
%   Copyright 2009-2010 Trustees of Columbia University

function h = cross_fb(t,n1,n2,tauf,scale)

h = scale/tauf*exp(-t/tauf).*((t/tauf).^n1/factorial(n1)-(t/tauf).^n2/factorial(n2));