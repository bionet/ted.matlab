%% variational example - single LIF neuron

clear all
close all

%% construct input

dt=1e-6;
I12=cumsum(ones(1,12));
Omega=2*pi*100;
T=pi/Omega;
t=0*T:dt:40*T;
Ik=[zeros(1,2), 2*(rand(1,36)-0.5), zeros(1,2)];

A=ones(length(Ik),1)*t-(1:length(Ik))'*T*ones(1,length(t));
        %A is representing (t-KT)
B=sinc(Omega/pi*A);
u=sum((Ik'*ones(1,length(t))).*B,1);

figure;plot(t,u); xlabel('time'); ylabel('current'); 
title('Bandlimited Input Stimuli');

%% single leaky neuron

b = 3;
d = 0.8;
C = 0.01;
R = 50;

tk  = iaf_encode_leaky(u, t, b, d, R, C);

u_rec = consistent_decoding_LIF(tk,t, b, d, R, C);

%% plot

figure;plot(t,u,t,u_rec); 
    xlabel('Time [sec]'); ylabel('Amplitude'); 
    title('Reconstruction of Temporal Contrast');
    legend('Original','Consistent Rec.');
    
t_v = round(0.1*length(t)):round(0.9*length(t));

snr = 10*log10(sum(u(t_v).^2)/sum((u(t_v)-u_rec(t_v)).^2)) % SNR for consistent recovery