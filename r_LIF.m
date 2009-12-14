function r = r_LIF(tk,RC)

r = (RC^2)*(tk(2:end)/RC-1-(tk(1:end-1)/RC-1).*exp(-diff(tk)/RC))';
