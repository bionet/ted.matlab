function tk  = iaf_encode_perfect(u, t, b, d, k)
% time encode with a perfect IAF neuron:
% u: the input signal
% t: time vector
%   and parameters:
% d: delta threshold of the IAF neuron
% b: bias 
% k: capacitance 
% and y(0) = 0

% and output
% tk: the spike times in seconds

% Time Encoding Tutorial Code
% BIONET Research Group, 09/2006

y = zeros(1, length(t)); %output of integrator (initialization)
s = zeros(1, length(t)); %marker for trigger times (initialization)
t0 = 1; 
dt = t(2)-t(1);
y(1)=dt*(b+u(1));
for i=2:length(t)    
    y(i)=y(i-1)+dt*(1/k)*(b+u(i));
    if y(i) >= d
        s(i) = 1;
        y(i) = y(i)-d; %reseting of the integrator to zero
        t0 = i;
    end
end

%calculation of spike times
tk=t(find(s));

%figure;  %uncheck to print
%subplot('position',[0.13 0.13 1-0.26 0.6])
%   plot(t,y);
%   hold on; plot([0 100],[10 10],'--');
%   axis([min(t) max(t) min(y) 1.1*max(y)])
%   xlabel('time [\tau]');
%   ylabel('u(t)')

%subplot('position',[0.13 0.8 1-0.26 0.1])
%   plot(t,s,'.','markersize',20);
%   axis([min(t) max(t) 0.5 1.5])
%   set(gca,'xtick',[],'ytick',[])
%   ylabel('spikes')
%title('Response of Time Encoding Machine');
