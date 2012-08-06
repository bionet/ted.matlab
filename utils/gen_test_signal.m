%GEN_TEST_SIGNAL Generate bandlimited, uniformly sampled test signal
%   X = GEN_TEST_SIGNAL(DUR,DT,FMAX) returns a signal DUR seconds
%   long uniformly sampled at intervals DT. The generated signal's 
%   frequency components are no higher than FMAX Hz.
%
%   X = GEN_TEST_SIGNAL(DUR,DT,FMAX,NP) adds Gaussian white noise
%   with power NP (in dB) to the generated signal before filtering
%   it. If NP is -inf, no noise will be added.
%
%   X = GEN_TEST_SIGNAL(DUR,DT,FMAX,NP,NC) sets the number of
%   frequency components in the generated signal to NC (the default is
%   3). If NC is 0, the generated signal will either be zero if NP
%   is -inf or pure filtered Gaussian white noise otherwise.

%   Author: Lev Givon
%   Copyright 2009-2012 Lev Givon

function u = gen_test_signal(dur,dt,fmax,varargin)

% The maximum frequency may not exceed the Nyquist frequency:
fs = 1/dt;
if fmax > fs/2,
  error('maximum frequency may not exceed the Nyquist frequency');
end

% Process optional parameters
np = -inf;
nc = 3;
if nargin > 3,
  np = varargin{1};
end
if nargin > 4,
  nc = varargin{2};
end

% Determine number of entries in generated signal. This corresponds
% to the length of [0:dt:dur]:
n = fix(dur/dt)+1;

% Randomly set nc frequency components:
f = zeros(1,n);
fmaxi = floor(n*fmax/fs);
if fmaxi < nc,
  error('maximum frequency is too low to provide %i frequency components',nc);
end

% The first element in the fft corresponds to the DC component; 
% hence, it is not set:
ci = [];
while length(ci) < nc,
  temp = randint(2,fmaxi+2);
  while intersect(ci,temp),
    temp = randint(2,fmaxi+2);
  end
  ci = [ci, temp];
end
p = -2*pi*rand(1,nc);
f(ci) = (n/2)*exp(1j*p);
%f(end-ci+2) = (n/2)*exp(1j*-p);

% Create the signal by transforming the constructed frequency 
% representation into the time domain and adding white noise if so
% specified:
u = real(ifft(f)) + randn(1,n)*10^(np/20);

% Filter the result. Since a cutoff of 1 corresponds to the Nyquist
% frequency 1/(2*dt), the cutoff corresponding to the frequency
% fmax must be fmax/(1/2*dt):
b = fir1(40,2*fmax*dt);
u = filter(b,1,u);

% Generate a random integer in the range [a,b):
function r = randint(a,b)

if a >= b,
  error('maximum bound must exceed minimum bound');
end
r = a + floor((b-a)*rand());

