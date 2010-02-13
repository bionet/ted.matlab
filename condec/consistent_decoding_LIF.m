function u_rec = consistent_decoding_LIF(tk,t, b, d, R, C)
% (tk,t,dt,b,d,R,C)
% consistent_decoding_LIF reconstruct signals encoded with a LIF neuron

% u_rec = consistent_decoding_LIF(t, b, d, R, C, lamda) reconstructs the encoded 
% stimulus u that has finite energy and is encoded with a 
% LIF neuron. The process is described in detail in section 2.2

% Inputs:
% tk:     spike times
% t:      time vector
% b:      bias
% delta:  threshold of the LIF neuron
% R:      resistance 
% C:      capacitance

% Output:
% u_rec:  reconstructed stimulus.

% Author(s): Eftychios A. Pnevmatikakis
% Copyright 2010 Trustees of Columbia University

ln = length(tk)-1;

RC=R*C;

q = C*(d-b*R)+ b*RC*exp(-(tk(2:end)-tk(1:end-1))/RC)';

G = G_IF(tk,RC);
p = p_IF(tk',RC);
r = r_IF(tk',RC);

V = [G p r; p' 0 0; r' 0 0];

cv = pinv(V)*[q;0;0];
d0 = cv(end-1);
d1 = cv(end);
c  = cv(1:end-2);

u_rec = d0+d1*t;

for i = 1:ln
    u_rec = u_rec + c(i)*psi(tk(i),tk(i+1),t,RC);
end