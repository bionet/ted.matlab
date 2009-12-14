function [tk1 tk2]  = iaf_encode_ideal_on_off(u, t, b, d, k, tauf, scale)
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

i=1;
y1 = zeros(1, length(t)); %output of integrator (initialization)
s1 = zeros(1, length(t)); %marker for trigger times (initialization)
v1 = u;

y2 = zeros(1, length(t)); %output of integrator (initialization)
s2 = zeros(1, length(t)); %marker for trigger times (initialization)
v2 = u;

t0 = 1; 
dt = t(2)-t(1);
y1(1)=dt*(b+v1(1));
y2(1)=dt*(-b+v2(i));

for i=2:length(t)    
    y1(i)=y1(i-1)+dt*(1/k)*(b+v1(i));
    y2(i)=y2(i-1)+dt*(1/k)*(-b+v2(i));
    
    if y1(i) >= d
        s1(i) = 1;
        y1(i) = y1(i)-d; %reseting of the integrator to zero
        %v1(i:end) = v1(i:end) + B*exp(-(t(i:end)-t(i))/tau);
        v2(i:end) = v2(i:end) - cross_fb(t(i:end)-t(i),5,7,tauf,scale); 
    end
    
    if y2(i) <= -d
        s2(i) = 1;
        y2(i) = y2(i)+d; %reseting of the integrator to zero
        %v2(i:end) = v2(i:end) - B*exp(-(t(i:end)-t(i))/tau);
        v1(i:end) = v1(i:end) + cross_fb(t(i:end)-t(i),5,7,tauf,scale);%[0,Bf*exp(-(t(i+1:end)-t(i+1))/tau_f)];        
    end    
    
end

%calculation of spike times
tk1 = t(find(s1));
tk2 = t(find(s2));

figure;plot(t,u,t,y1,t,y2); hold on;
xlabel('Time [sec]'); ylabel('Amplitude');
plot(tk1-dt,d*ones(1,length(tk1)),'yo');hold on;%y1(find(s1)),'yo'); hold on
plot(tk2-dt,-d*ones(1,length(tk2)),'mo');hold on;%y2(find(s2)),'mo'); hold on
plot(t,d*ones(size(t)),'k--',t,-d*ones(size(t)),'k--'); hold off;
legend('Input','v^1','v^2','ON spikes','OFF spikes');
title(['ON-OFF IAF \delta= ' num2str(d)]);