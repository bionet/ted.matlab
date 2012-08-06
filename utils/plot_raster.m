%PLOT_RASTER Create raster plot of time sequences.
%   PLOT_RASTER({T1,T2,...}) plots the specified time sequences as
%   rasters. The sequences must be specified as elements of a cell array.

%   Author: Lev Givon
%   Copyright 2009-2012 Lev Givon

function plot_raster(t_list)

clf;
x_max = 0;
n = length(t_list);
for i=1:n,
    if max(t_list{i}) > x_max,
        x_max = max(t_list{i});
    end
end
y_offset = 0.5;
ha = axes();
axis(ha,[0,x_max,y_offset,y_offset+n]);
set(ha,'ytick',[1:n]);
set(ha,'yticklabel',[1:n]);
set(ha,'box','on');
gap = 0.1;
for i=1:length(t_list),
    t = t_list{i};
    for j=1:length(t),
        hl = line([t(j),t(j)],[y_offset+gap,y_offset+1-gap]);
        set(hl,'linewidth',1);
    end
    hl = line([0,x_max],[y_offset+1,y_offset+1]);
    set(hl,'color','k');
    y_offset = y_offset+1;
end
xlabel('t (seconds)');
