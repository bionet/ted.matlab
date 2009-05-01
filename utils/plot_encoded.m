%PLOT_ENCODED Plot results of time encoding.
%   PLOT_ENCODED(T,U,S) plots the time-encoded spike train S
%   alongside the original signal U over the uniformly spaced times
%   T. S is assumed to contain the intervals between spikes.
%
%   PLOT_ENCODED(T,U,S,FIG_TITLE) plots U and S over T with the
%   plot title FIG_TITLE.
%
%   PLOT_ENCODED(T,U,S,FIG_TITLE,FILENAME) saves the plot in
%   FILENAME.

function plot_encoded(t,u,s,varargin)

if nargin > 5,
  error('Too many input arguments.');
end

clf;
dt = t(2)-t(1);
cs = cumsum(s);
a = axes('position',[0.13,0.3,0.775,0.6150]);
stem(cs,u(round(cs/dt)));
hold on;
plot(t,u);
xlabel('t (seconds)');
ylabel('u(t)');

if nargin >= 4,
  title(varargin{1});
end

set(a,'XLim',[min(t),max(t)]);
a = axes('position',[0.13,0.1,0.775,0.1]);
plot(cs,zeros(1,length(s)),'ro');
xlabel(sprintf('%d spikes',length(s)));
set(a,'XLim',[min(t),max(t)]);
set(a,'YTick',[]);

if nargin == 5,
  filename = varargin{2};
  tok = regexp(filename,'\.(?<first>.*)$','tokens');
  if isempty(tok),
    error('Unrecognized output file format.');
  end
  filetype = tok{1}{1};
  print(gcf,strcat('-d',filetype),filename);
end
