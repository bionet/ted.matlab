%IAF_ENCODE Encode a signal using an integrate-and-fire neuron
%   S = IAF_ENCODE(U,DT,B,D,SD,R,C,DTE,Y,INT,QM,FOUT) encodes the
%   signal U sampled at rate 1/DT using an integrate-and-fire neuron
%   with bias B, Gaussian distributed firing threshold with mean D, 
%   standard deviation SD, resistance R, and capacitance C. If
%   DTE is specified and is smaller than DT, U is resampled with
%   at the sampling frequency 1/DTE Hz prior to being integrated. The
%   state of the IAF integrator is initially set to Y and the initial
%   interval since the last spike is set to INT. One may specify the
%   quadrature method used to integrate U; currently 'rect'
%   (rectangular integration) and 'trapz' (trapezoidal integration)
%   are supported. The encoded signal is returned as an array of
%   intervals S. 
%
%   The parameters K, DTE, Y, INT, QM, and FOUT are optional. If R =
%   inf (the default), an ideal neuron model is used. C is assumed
%   to be 1 if not specified.
%
%   [S,Y,INT] = IAF_ENCODE(U,T,B,D,R,C,DTE,Y,INT,QM,FOUT) returns
%   the values of the state Y of the IAF neuron, and the time INT since
%   the last emitted spike.

%   Author: Lev Givon
%   Copyright 2009-2014 Lev Givon

function [varargout] = iaf_encode(u,dt,b,d,varargin)

R = inf;
C = 1;
dte = dt;
y = 0;
interval = 0;
quad_method = 'trapz';
full_output = false;

if nargin >= 5, 
  sd = varargin{1};
else
  sd = 0;
end
if nargin >= 6, 
  R = varargin{2};
end
if nargin >= 7,
  C = varargin{3};
end
if nargin >= 8,
  dte = varargin{4};
end
if nargin >= 9,
  y = varargin{5};
end
if nargin >= 10,
  interval = varargin{6};
end
if nargin >= 11,
  quad_method = varargin{7};
end
if nargin >= 12,
  full_output = varargin{8};
end
if nargin >= 13,
  error('Too many input arguments.');
end

if isempty(u),
  if full_output,
    varargout = {[],y,interval};
  else
    varargout = {[]};
  end
  return
end

% Check whether the encoding resolution is finer than that of the
% original sampled signal:
nu = length(u);
if dte > dt,
  error(['Encoding time resolution must not exceed original signal ' ...
         'resolution.']);
end
if dte < 0,
  error('Encoding time resolution must be nonnegative.');
end
if dte ~= 0 && dte ~= dt,
  u = resample(u,round(dt/dte),1);
  nu = length(u);
  dt = dte;
end

s = [];

% Choose integration method:
if isinf(R),
  if strcmp(quad_method,'rect'),
    compute_y = @(y,i) y + dt*(b+u(i))/C;
    last = nu;
  elseif strcmp(quad_method,'trapz'),
    compute_y = @(y,i) y + dt*(b+(u(i)+u(i+1))/2)/C;
    last = nu-1;
  else
    error('Unrecognized quadrature method.');
  end
else
  
  % When the neuron is leaky, use the exponential Euler method to
  % perform the encoding:
  RC = R*C;
  compute_y = @(y,i) y*exp(-dt/RC)+R*(1-exp(-dt/RC))*(b+u(i));
  last = nu;
end

% The interval between spikes is saved between iterations rather than
% the absolute time so as to avoid overflow problems for very long
% signals:
j = 1;

d_new = d + sd * randn;
while d_new < 0  % force positive threshold
  d_new = d + sd * randn;
end

for i=1:last,
  y = compute_y(y,i);
  interval = interval + dt;
  if y >= d_new
    s(j) = interval;
    j = j + 1;
    interval = 0;
    y = y - d_new;
    d_new = d + sd * randn;
    while d_new < 0
      d_new = d + sd * randn;
    end
  end
end

if full_output,
  varargout = {s,y,interval};
else
  varargout = {s};
end
