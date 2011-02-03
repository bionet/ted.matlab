%CONSISTENT_DECODING_LIF Decode a signal encoded with a leaky IAF neuron.
%   U_REC = CONSISTENT_DECODING_LIF(TK,T,B,D,R,C) decodes a
%   finite-energy signal encoded as spike times TK over the times T
%   using a leaky IAF neuron with bias B, threshold D, resistance R,
%   and capacitance C.
%
%   The calculation is described in further detail in Section 2.2 of
%   the Consistent Recovery paper mentioned in the toolbox references.

%   Author: Eftychios A. Pnevmatikakis
%   Copyright 2009-2011 Eftychios A. Pnevmatikakis

function u_rec = consistent_decoding_LIF(tk,t,b,d,R,C)

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
