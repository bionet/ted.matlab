function psi = psi_LIF(tk1,tk2,RC,t)

psi = zeros(1,length(t));

t1 = (tk1 - t)/RC;
t2 = (tk2 - t)/RC;

f=@(x)(x.^3-3*x.^2+6*x-6);

ft  = (RC^4)*f([t1;t2]);
etk = exp(-(tk2-tk1)/RC);
ett = 12*(RC^4)*exp(-[t1;t2]);

fp = find(t>tk2,1);
fm = find(t>tk1,1);

psi(1:fm-1)  = ft(2,1:fm-1)-ft(1,1:fm-1)*etk;
psi(fm:fp-1) = ett(2,fm:fp-1)+ft(2,fm:fp-1)+ft(1,fm:fp-1)*etk;
psi(fp:end)  = ft(1,fp:end)*etk-ft(2,fp:end);