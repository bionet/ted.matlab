%FFTCONV Convolve two signals using FFT.
%   Y = FFTCONV(A,B) convolves vectors A and B. It uses the FFT to
%   accelerate the computation.

%   Authors: Eftychios A. Pnevmatikakis and Robert J. Turetsky
%   Copyright 2009-2011 Lev Givon

function y = fftconv(a,b)

y = real(ifft(fft([a zeros(1,length(b)-1)]).* ...
              fft([b zeros(1,length(a)-1)])));

