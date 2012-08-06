%ASDM_ENCODE Encode a signal with an asynchronous sigma-delta modulator.
%   S = ASDM_ENCODE(U,DT,B,D,K,DTE,Y,INT,SGN,QM,FOUT) encodes the
%   signal U sampled at rate 1/DT using an asynchronous sigma-delta
%   modulator with bias B, Schmitt trigger threshold D, and
%   integration constant K. If DTE is specified and is smaller than
%   DT, U is resampled with at the sampling frequency 1/DTE Hz prior
%   to being integrated. The state of the ASDM integrator is initially
%   set to Y, the initial interval since the last spike is set to INT,
%   and the sign of the integrator at the first spike is set to
%   SGN. One may specify the quadrature method used to integrate U;
%   currently 'rect' (rectangular integration) and 'trapz'
%   (trapezoidal integration) are supported. The encoded signal is
%   returned as an array of intervals S. K, DTE, Y, INT, QM, and FOUT
%   are optional.
%
%   [S,Y,INT,SGN] = ASDM_ENCODE(U,DT,B,D,K,DTE,Y,INT,QM,FOUT)
%   returns the values of the state Y of the ASDM integrator, the time
%   INT since the last emitted spike, and the sign SGN of the
%   integrator at the conclusion of processing U.

%   Author: Lev Givon
%   Copyright 2009-2012 Lev Givon

function [varargout] = asdm_encode(u,dt,b,d,varargin)

k = 1;
dte = dt;
y = 0;
interval = 0;
sgn = 1;
quad_method = 'trapz';
full_output = false;

if nargin >= 5,
  k = varargin{1};
end
if nargin >= 6,
  dte = varargin{2};
end
if nargin >= 7,
  y = varargin{3};
end
if nargin >= 8,
  interval = varargin{4};
end
if nargin >= 9,
  sgn = varargin{5};
end
if nargin >= 10,
  quad_method = varargin{6};
end
if nargin >= 11,
  full_output = varargin{7};
end
if nargin >= 12,
  error('Too many input arguments.');
end

if isempty(u),
  if full_output,
    varargout = {[],y,interval,sgn};
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
if (dte ~= 0) && (dte ~= dt),
  u = resample(u,round(dt/dte),1);
  nu = length(u);
  dt = dte;
end

s = [];

% Choose integration method:
if strcmp(quad_method,'rect'),
  compute_y = @(y,sgn,i) y + dt*(sgn*b+u(i))/k;
  last = nu;
elseif strcmp(quad_method,'trapz'),
  compute_y = @(y,sgn,i) y + dt*(sgn*b+(u(i)+u(i+1))/2)/k;
  last = nu - 1;
else
  error('Unrecognized quadrature method.');
end

% The interval between spikes is saved between iterations rather than
% the absolute time so as to avoid overflow problems for very long
% signals:
j = 1;
for i=1:last,
  y = compute_y(y,sgn,i);
  interval = interval + dt;
  if abs(y)>=d,
    s(j) = interval;
    j = j + 1;
    interval = 0;
    y = d*sgn;
    sgn = -sgn;
  end
end

if full_output,
  varargout = {s,y,interval,sgn};
else
  varargout = {s};
end
