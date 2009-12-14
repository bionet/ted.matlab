function u_rec = consistent_decoding_IF_ONOFF(tk1, tk2, t, b, d, C, tauf, scale)

dt = t(2) - t(1);
ln = [length(tk1)-1,length(tk2)-1];
ln2 = cumsum([0,ln]);
q = q_IF_ONOFF(tk1,tk2, b, d, C, tauf, scale, dt);

TK(1:ln(1)+1,1)=tk1';
TK(1:ln(2)+1,2)=tk2';

G = G_pop_IF(TK,ln,2);
p = p_pop_IF(TK,ln,2);
r = r_pop_IF(TK,ln,2);

V = [G p r; p' 0 0; r' 0 0];

cv = pinv(V)*[q;0;0];
d0 = cv(end-1);
d1 = cv(end);
c  = cv(1:end-2);

u_rec = d0+d1*t;

for j = 1:2
    for i = 1:ln(j)
        u_rec = u_rec + c(ln2(j)+i)*psi_IF(TK(i,j),TK(i+1,j),t);
    end
end