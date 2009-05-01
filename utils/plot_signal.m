%PLOT_SIGNAL Plot a signal.
%   PLOT_SIGNAL(T,U) plots the signal U over the
%   uniformly spaced times T. 
%
%   PLOT_SIGNAL(T,U,FIG_TITLE) plots U over T with the plot title 
%   FIG_TITLE. 
%  
%   PLOT_SIGNAL(T,U,FIG_TITLE,FILENAME) saves the plot in FILENAME.

function plot_signal(t,u,varargin)

if nargin > 4,
  error('Too many input arguments.');
end

clf;
plot(t,u);
xlabel('t (seconds)'); 
ylabel('u(t)');
set(gca,'XLim',[min(t),max(t)]);

if nargin >= 3,
  set(gcf,'Name',varargin{1});
  title(varargin{1});
end

if nargin == 4,
  filename = varargin{2};
  tok = regexp(filename,'\.(?<first>.*)$','tokens');
  if isempty(tok),
    error('Unrecognized output file format.');
  end
  filetype = tok{1}{1};
  print(gcf,strcat('-d',filetype),filename);
end

