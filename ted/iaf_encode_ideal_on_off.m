%IAF_ENCODE_IDEAL_ON_OFF Encode a signal with an ON-OFF IAF neuron pair.
%   [TK1 TK2] = IAF_ENCODE_IDEAL_ON_OFF(U,T,B,D,K,TAUF,SCALE) encodes
%   a signal U defined over the times T using an ON-OFF IAF neuron
%   pair with bias B, threshold D and integration constant K. The
%   neurons are coupled with a cross-feedback term described by the
%   time constant TAUF and the amplitude scaling factor SCALE.
%
%   The calculation is described in further detail in Equation 33
%   of the Consistent Recovery paper mentioned in the toolbox references.

%   Author: Eftychios A. Pnevmatikakis
%   Copyright 2009-2010 Trustees of Columbia University

function [tk1 tk2]  = iaf_encode_ideal_on_off(u,t,b,d,k,tauf,scale)

i = 1;
y1 = zeros(1,length(t)); % initialize output of integrator 
s1 = zeros(1,length(t)); % initialize trigger time marker
v1 = u;

y2 = zeros(1,length(t)); % initialize output of integrator 
s2 = zeros(1,length(t)); % initialize trigger time marker
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
        y1(i) = y1(i)-d; % reset the integrator to 0
        % v1(i:end) = v1(i:end) + B*exp(-(t(i:end)-t(i))/tau);
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
xlabel('t (seconds)'); ylabel('u(t)');
plot(tk1-dt,d*ones(1,length(tk1)),'yo');hold on;%y1(find(s1)),'yo'); hold on
plot(tk2-dt,-d*ones(1,length(tk2)),'mo');hold on;%y2(find(s2)),'mo'); hold on
plot(t,d*ones(size(t)),'k--',t,-d*ones(size(t)),'k--'); hold off;
legend('Input','v^1','v^2','ON spikes','OFF spikes');
title(['ON-OFF IAF \delta= ' num2str(d)]);