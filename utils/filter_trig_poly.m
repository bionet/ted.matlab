%FILTER_TRIG_POLY Filter a trigonometric polynomial
%   Y = FILTER_TRIG_POLY(U,H) filters the trigonometric polynomial
%   signal U with a filter whose impulse response is given in H
%   using the FFT.

%   Author: Lev Givon 
%   Copyright 2009-2012 Lev Givon

function y = filter_trig_poly(u,h)

Nu = length(u);
Nh = length(h);
y = real(ifft(fft([u zeros(1,Nh-1)]).*fft([h zeros(1,Nu-1)])));
y = y(1:Nu);
