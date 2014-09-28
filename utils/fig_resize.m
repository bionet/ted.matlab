%FIG_RESIZE Resize a figure.
%   FIG_RESIZE(H,R) resizes the figure specified by handle H by
%   rescaling the length and width by R.
%
%   FIG_RESIZE(H,RW,RL) resizes the figure's width by RW and length
%   by RL.

%   Author: Lev Givon
%   Copyright 2009-2014 Lev Givon

function fig_resize(h,varargin)

if nargin > 3,
    error('Too many arguments.');
end

rw = varargin{1};
rl = rw;
if nargin == 3,
    rl = varargin{2};
end

pos = get(h,'Position');
pos(3) = rw*pos(3);
pos(4) = rl*pos(4);
set(h,'Position',pos);
