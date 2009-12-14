function tk  = iaf_encode_leaky(u, t, b, d, R, C)
% time encode with leaky IAF neuron:
% x: the input signal
% t: time vector
%   and parameters:
% d: threshold of the IAF neuron
% b: bias
% R: resistance 
% C: capacitance
% and output
% tk: the spike times in seconds

% Time Encoding Tutorial Code
% BIONET Research Group, 09/2006

y = zeros(1, length(t)); %output of integrator (initialization)
s = zeros(1, length(t)); %marker for trigger times (initialization)
t0 = 1; 
dt = t(2)-t(1);
y(1)=dt*(b+u(1))*exp(-dt/R/C);
for i=2:length(t)
    y(i) = y(i-1)*exp(-dt/R/C)+R*(1-exp(-dt/R/C))*(u(i)+b);
    if y(i) >= d
        s(i) = 1;
        y(i) = y(i)-d; %reseting of the integrator to zero
        t0 = i;
    end
end

%calculation of spike times
tk=t(find(s));