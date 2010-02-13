% IF_ONOFF_demo
% Performs the example presented in section 3.3 
% The input signal represents the temporal contrast of a positive input
% photocurrent. It is encoded with a pair of interconnected ON-OFF IF neuron. 
% The signal is then reconstructed within the L2 space, without assuming 
% any prior knowledge

% Author: Eftychios A. Pnevmatikakis
% Copyright 2010 Trustees of Columbia University

% run time ~ 40sec on an iMac 5.1, Intel Core 2 Duo, 2.16GHz, 4GB RAM

%% construct temporal contrast input

dt = 1e-6;
fmax = 100;
t = dt:dt:0.2;
tr_vc = round(0.05*length(t)):round(0.95*length(t)); % truncated vector

b_p = 8;

mc = floor(floor(t(end)/dt)*fmax*dt); % maximum sinusoidal components
v = b_p + gen_test_signal(t(end)+2*round(0.1*length(t))*dt,dt,fmax,-Inf,mc);
    % a constant bias is added to the signal to ensure that it is positive
    
v = v((round(0.1*length(t))+1):end-round(0.1*length(t))); % truncate first 100 values to eliminate discontinuities   


figure;subplot(1,2,1);plot([0,t],v)
    xlabel('Time [sec]');
    ylabel('Intensity');
    title('Input Stimulus I(t)')

u = (diff(v)/dt)./v(2:end);
u = u/max(abs(u));

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

figure;plot(t(tr_vc),u(tr_vc),t(tr_vc),u_rec(tr_vc)); 
    xlabel('Time [sec]'); ylabel('Amplitude'); 
    title('Reconstruction of Temporal Contrast');
    legend('Original','Consistent Rec.');
    
snr = 10*log10(sum(u(tr_vc).^2)/sum((u(tr_vc)-u_rec(tr_vc)).^2)) % SNR for consistent recovery