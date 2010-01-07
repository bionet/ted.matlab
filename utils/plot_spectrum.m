%PLOT_SPECTRUM Plot a signal's frequency spectrum.
%   PLOT_SPECTRUM(U,FS) plots the frequency spectrum of the signal
%   U assuming a sampling frequency of FS Hz. 
%
%   PLOT_SPECTRUM(U,FS,FMIN) plots the portion of the spectrum
%   between FMIN and FS/2 Hz.
%
%   PLOT_SPECTRUM(U,FS,FMIN,FMAX) plots the portion of the spectrum
%   between FMIN and FMAX Hz.

%   Author: Lev Givon
%   Copyright 2009-2010 Trustees of Columbia University

function plot_spectrum(u,fs,varargin)

if nargin > 2,
  fmin = varargin{1};
else
  fmin = 0;
end
if fmin < 0 || fmin >= fs/2,
  error('invalid minimum frequency');
end

if nargin == 4,
  fmax = varargin{2};
else
  fmax = fs/2;
end
if fmax <= fmin || fmax > fs/2,
  error('invalid maximum frequency');
end

n = floor(length(u)/2);
uf = fft(u);
uf = uf(1:n);
f = (fs/2)*[0:n-1]/n;

a = floor(2*n*fmin/fs)+1;
b = floor(2*n*fmax/fs);

length(uf)
clf();
subplot(211);
stem(f(a:b),real(uf(a:b)));
ylabel('real');
subplot(212);
stem(f(a:b),imag(uf(a:b)));
ylabel('imag');
xlabel('f (Hz)');
