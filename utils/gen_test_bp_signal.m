%GEN_TEST_BP_SIGNAL Generate bandlimited, uniformly sampled test signal
%   X = GEN_TEST_BP_SIGNAL(DUR,DT,FMIN,FMAX) returns a signal DUR
%   seconds long uniformly sampled at intervals DT. The generated
%   signal's frequency components are no lower than FMIN Hz and no
%   higher than FMAX Hz.
%
%   X = GEN_TEST_BP_SIGNAL(DUR,DT,FMIN,FMAX,NP) adds Gaussian white
%   noise with power NP (in dB) to the generated signal before
%   filtering it. If NP is -inf, no noise will be added.
%
%   X = GEN_TEST_BP_SIGNAL(DUR,DT,FMIN,FMAX,NP,NC) sets the number of
%   frequency components in the generated signal to NC (the default is
%   3). If NC is 0, the generated signal will either be zero if NP is
%   -inf or filtered pure Gaussian white noise otherwise.

%   Authors: Lev Givon and Eftychios A. Pnevmatikakis 
%   Copyright 2009-2015 Lev Givon and Eftychios A. Pnevmatikakis

function u = gen_test_bp_signal(dur,dt,fmin,fmax,varargin)

% The minimum frequency may not be negative:
if fmin < 0
    error('minimum frequency may not be negative');
end

% The minimum frequency may not exceed the maximum frequency:
if fmin > fmax
    error('minimum frequency may not exceed maximum frequency');
end

% The maximum frequency may not exceed the Nyquist frequency:
fs = 1/dt;
if fmax > fs/2,
  error('maximum frequency may not exceed the Nyquist frequency');
end

% Process optional parameters
np = -inf;
nc = 3;
if nargin > 4,
  np = varargin{1};
end
if nargin > 5,
  nc = varargin{2};
end

n = round(dur/dt);

% Randomly set nc frequency components:
f = zeros(1,n);
fmaxi = floor(n*fmax/fs);
fmini = ceil(n*fmin/fs);
if fmaxi-fmini-1 < nc,
  error('frequency support is too narrow to provide %i frequency components',nc);
end

ci = [];
while length(ci) < nc,
  temp = randi([fmini+2,fmaxi+1]);
  while intersect(ci,temp),
    temp = randi([fmini+2,fmaxi+1]);
  end
  ci = [ci, temp];
end
p = -pi*rand(1,nc);
f(ci) = (n/2)*exp(1j*p);

% Create the signal by transforming the constructed frequency 
% representation into the time domain and adding white noise if so
% specified:
u = real(ifft(f))+randn(1,n)*10^(np/20);

% Filter the result. Since a cutoff of 1 corresponds to the Nyquist
% frequency 1/(2*dt), the cutoff corresponding to the frequency
% fmax must be fmax/(1/2*dt):
b = fir1(40,[2*fmin*dt,2*fmax*dt]);
u = filter(b,1,u);
u = u(1:n);
