function h = cross_fb(t,n1,n2,tauf,scale)

h = scale/tauf*exp(-t/tauf).*((t/tauf).^n1/factorial(n1)-(t/tauf).^n2/factorial(n2));