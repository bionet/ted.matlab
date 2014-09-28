%MAKE_GAMMATONE_FB Gammatone filterbank generator
%   H = MAKE_GAMMATONE_FB(T,N) computes the impulse responses H
%   of a gammatone filterbank of N filters. Each filter's impulse
%   response  is stored in a row of H. The filters are all
%   normalized such that the maximum magnitude of each filter's
%   spectrum is 0 dB.

%   Author: Lev Givon
%   Copyright 2009-2014 Lev Givon

function h = make_gammatone_fb(t,N)

h = zeros(N,length(t));
for i=1:N,
    h(i,:) = make_gammatone(t,N,i);
end

% Compute a filter in the filter bank:
function h = make_gammatone(t,N,i)

if i<1 || i>N,
    error('invalid value for i');
end

Q_ear = 9.26449;
BW_min = 24.7;
beta = 1.019;
n = 4;
f_min = 200;
f_max = 1400;

ERB = @(f_c) ((f_c/Q_ear)^n+(BW_min)^n)^(1/n);

o = (Q_ear/N)*(log(f_max+Q_ear*BW_min)-log(f_min+Q_ear*BW_min));
f_c = -Q_ear*BW_min+(f_max+Q_ear*BW_min)*exp(-(o/Q_ear)*(N-i));
h = (t.^(n-1)).*exp(-2*pi*beta*ERB(f_c).*t).*cos(2*pi*f_c*t);

% Scale the impulse response so that the magnitude of the largest 
% frequency components is 1:
h = h/max(abs(fft(h)));

