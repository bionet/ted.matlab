function p = p_LIF(tk,RC)

p = RC*(1-exp(-diff(tk)/RC))';
