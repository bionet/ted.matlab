%GEN_WINDOW Generate a tapered window
%   Y = GEN_WINDOW(N,O,Z) generates a tapered window with a middle
%   section comprising N points set to 1, tapered portions obtained
%   from a Hann window each comprising O points, and zero portions
%   each comprising Z points. The total length of the generated
%   window is N+2*O+2*Z. This is equivalent to invoking the
%   function as GEN_WINDOW(N,O,Z,'f').
%
%   Y = GEN_WINDOW(N,O,Z,'r') generates a window with a right-hand
%   tapered and zeroed portion constructed as above; the remaining
%   part of the window is set to 1. As before, the total length of
%   the generated window is N+2*O+2*Z.
%
%   Y = GEN_WINDOW(N,O,Z,'l') generates a window with a left-hand
%   tapered portion.

%   Author: Lev Givon
%   Copyright 2009-2014 Lev Givon

function y = gen_window(N,O,Z,varargin)

wintype = 'f';
if nargin == 4,
    wintype = varargin{1};
end

h = hann(O*2)';
if wintype == 'f',
    y = [zeros(1,Z), h(1:O), ones(1,N), h(O+1:end), zeros(1,Z)];
elseif wintype == 'r',
    y = [ones(1,N+Z+O), h(O+1:end), zeros(1,Z)];
elseif wintype == 'l',
    y = [zeros(1,Z), h(1:O), ones(1,N+Z+O)];
else
    error('invalid window type');
end
