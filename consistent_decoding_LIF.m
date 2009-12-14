function u_rec = consistent_decoding_LIF(tk,t, b, d, R, C)

ln = length(tk)-1;

RC=R*C;

q = C*(d-b*R)+ b*RC*exp(-(tk(2:end)-tk(1:end-1))/RC)';

G = G_LIF(tk,RC);
p = p_LIF(tk,RC);
r = r_LIF(tk,RC);

V = [G p r; p' 0 0; r' 0 0];

cv = pinv(V)*[q;0;0];
d0 = cv(end-1);
d1 = cv(end);
c  = cv(1:end-2);

u_rec = d0+d1*t;

for i = 1:ln
    u_rec = u_rec + c(i)*psi_LIF(tk(i),tk(i+1),RC,t);
end