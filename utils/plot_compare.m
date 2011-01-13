%PLOT_COMPARE Compare two signals and plot their differences.
%   PLOT_COMPARE(T,U,V) plots the signal U and V on the same
%   axis. It also plots the difference between them.

%   Author: Lev Givon
%   Copyright 2009-2011 Lev Givon

function plot_compare(t,u,v,varargin)

if nargin > 5,
  error('Too many arguments.');
end

clf;
subplot(211);
plot(t,u,'b',t,v,'r');
xlabel('t (seconds)'); 
ylabel('u(t)');
set(gca,'XLim',[min(t),max(t)]);
subplot(212);
plot(t,20*log10(abs(u-v)));
xlabel('t (seconds)'); ylabel('error (dB)');
set(gca,'XLim',[min(t),max(t)]);

if nargin >= 4,
  subplot(211);
  set(gcf,'Name',varargin{1});
  title(varargin{1});
end

if nargin >= 5,
  filename = varargin{2};
  tok = regexp(filename,'\.(?<first>.*)$','tokens');
  if isempty(tok),
    error('Unrecognized output file format.');
  end
  filetype = tok{1}{1};
  print(gcf,strcat('-d',filetype),filename);
end
