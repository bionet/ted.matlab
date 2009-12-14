% consistent decoding integrate-and-fire on-off example

clear all
close all

%% construct temporal contrast input

dt=1e-6;
Omega=2*pi*100;
T=pi/Omega;
t=-0*T:dt:40*T;
rand('seed',63) % same input signal

Ik=[zeros(1,4), exprnd(5,1,32), zeros(1,4)]; % create positive bandlimited input

A=ones(length(Ik),1)*t-(1:length(Ik))'*T*ones(1,length(t));
        %A is representing (t-KT)
B=sinc(Omega/pi*A);

clear A;

b1 = 4; % bias to drive signal away from zero
v= b1 + sum((Ik'*ones(1,length(t))).*B,1);
clear B;


figure;subplot(1,2,1);plot(t,v)
    xlabel('Time [sec]');
    ylabel('Intensity');
    title('Input Stimulus I(t)')

u = (diff(v)/dt)./v(2:end);
u = u/max(abs(u));

t = t(1:end-1);
subplot(1,2,2); plot(t,u)
    xlabel('Time [sec]');
    ylabel('Contrast');
    title('Temporal Contrast v(t)')
    

%% encode

b = 4;          % neuron bias
d = 0.75;       % neuron threshold
C = 0.01;       % integration constant
tauf = 0.015;   % a - constant for the feedback
scale = 0.005;  % s - scaling constant

[tk1 tk2]  = iaf_encode_ideal_on_off(u, t, b, d, C, tauf, scale);

%% decoding

u_rec = consistent_decoding_IF_ONOFF(tk1, tk2, t, b, d, C, tauf, scale);

%% plot

figure;plot(t,u,t,u_rec); 
    xlabel('Time [sec]'); ylabel('Amplitude'); 
    title('Reconstruction of Temporal Contrast');
    legend('Original','Consistent Rec.');
    
t_v = round(0.1*length(t)):round(0.9*length(t));

snr = 10*log10(sum(u(t_v).^2)/sum((u(t_v)-u_rec(t_v)).^2)) % SNR for consistent recovery