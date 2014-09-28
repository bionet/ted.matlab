%GAMMATONE Create a gammatone filter bank.
%   [H,FC,T,F] = GAMMATONE(NUM,LEN,FMIN,FMAX,FS) creates
%   a gammatone filter bank containing NUM filters of length LEN
%   samples. The center frequency of the first filter is FMIN Hz
%   and the maximum frequency of the filter bank is FMAX Hz.
%   The generated filter bank is returned as a matrix of size NUM x
%   LEN transfer functions, a vector FC of center frequencies, and
%   the times T and frequencies F over which the filters are
%   constructed. 
%
%   [H,FC,T,F] = GAMMATONE(NUM,LEN,FMIN,FMAX,FS,PAD_BW) pads the
%   highest band in the filter bank is by PAD_BW to avoid aliasing; by
%   default, PAD_BW == 2; if PAD_BW == 1, then FC(end) == FMAX.

%   Authors: Eftychios A. Pnevmatikakis and Robert J. Turetsky
%   Copyright 2009-2014 Eftychios A. Pnevmatikakis and Robert J. Turetsky

function [h,fc,t,f] = gammatone(num,len,fmin,fmax,fs,pad_bw)

if exist('pad_bw') ~= 1
    pad_bw = 2;
end

EarQ = 9.26449;
minBW = 24.7;
order = 4;
dt = 1/fs;
t = dt*(0:len-1);
f = (0:length(t)-1)/length(t)*fs;
beta = 1.019;

Wp = fmax;
fmax = EarQ*(Wp-pad_bw*beta*minBW)/(EarQ+pad_bw*beta);

overlap = EarQ*(log(fmax+EarQ*minBW)-log(fmin+EarQ*minBW))/max(1,num-1);
fc = -EarQ*minBW + (fmax+EarQ*minBW)*exp(-(num-(1:num))*overlap/EarQ);
h = zeros(num,len);
for i=1:num,
    h(i,:) = t.^(order-1).*exp(-2*pi*beta*(fc(i)/EarQ+minBW)*t).*cos(2*pi*fc(i)*t);
    h(i,:) = h(i,:)/max(abs(fft(h(i,:))));
end
