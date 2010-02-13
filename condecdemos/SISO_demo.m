% SISO_demo
% Performs the example presented in section 2.3 
% The input signal is bandlimited with bandwidth 100Hz, restristed to 200ms
% It is encoded with a single LIF neuron. The signal is then reconstructed 
% within the L2 space, without taking advantege of the known bandwidth

% Author: Eftychios A. Pnevmatikakis
% Copyright 2010 Trustees of Columbia University

% run time ~ 8sec on an iMac 5.1, Intel Core 2 Duo, 2.16GHz, 4GB RAM

%% construct stimulus
dt = 1e-6;
fmax = 100;
t = dt:dt:0.2;
tr_vc = round(0.05*length(t)):round(0.95*length(t)); % truncated vector

mc = floor(floor(t(end)/dt)*fmax*dt); % maximum sinusoidal components
u = gen_test_signal(t(end)+(2*round(0.1*length(t))-1)*dt,dt,fmax,-Inf,mc);
    % a constant bias is added to the signal to ensure that it is positive
    
u = u((round(0.1*length(t))+1):end-round(0.1*length(t))); % truncate first 100 values to eliminate discontinuities 

u = 1.5*u/max(abs(u)); % normalization

%figure;plot(t,u); xlabel('time'); ylabel('current'); 
%title('Bandlimited Input Stimuli');

%% single leaky neuron

b = 3;
d = 0.8;
C = 0.01;
R = 50;

tk  = iaf_encode(u,dt,b,d,R,C);

u_rec = consistent_decoding_LIF(cumsum([0,tk]),t, b, d, R, C);
%% plot

figure;plot(t(tr_vc),u(tr_vc),t(tr_vc),u_rec(tr_vc)); 
    xlabel('Time [sec]'); ylabel('Amplitude'); 
    title('Reconstruction of Bandlimited Input');
    legend('Original','Consistent Rec.');
    
snr = 10*log10(sum(u(tr_vc).^2)/sum((u(tr_vc)-u_rec(tr_vc)).^2)) % SNR for consistent recovery