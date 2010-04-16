%% Time Encoding and Decoding Using ON/OFF Integrate-and-Fire Neurons
% This demo illustrates the time encoding of a bandlimited signal
% using interconnected ON/OFF integrate-and-fire neurons and recovery
% of the signal using spline interpolation. The input signal in this
% case can be thought of as the temporal contrast of a positive input
% photocurrent.
%
% The demo corresponds to the example presented in Section 3.3 of
% the Consistent Recovery paper mentioned in the toolbox
% references.

%% Generating a Test Signal
% Generate a noiseless signal 0.2 s long sampled at 1 GHz with a
% bandwidth of 100 Hz:
dur = 0.2;       % duration
dt = 1e-6;       % sampling resolution
fmax = 100;      % bandwidth (Hz)
t = [dt:dt:dur]; % time support

%%
% Truncate the time vector:
tr_vc = round(0.05*length(t)):round(0.95*length(t));

%%
% Add a constant bias to the signal to ensure that it is positive:
b_p = 8;

%%
% Use the maximum possible number of frequency components:
mc = floor(floor(t(end)/dt)*fmax*dt); 
rand('twister',0); randn('state',0);
v = b_p + gen_test_signal(t(end)+2*round(0.1*length(t))*dt,dt,fmax,-Inf,mc);
    
%%    
% Truncate first 100 values to eliminate discontinuities:
v = v((round(0.1*length(t))+1):end-round(0.1*length(t))); 

%%
% Normalize the generated signal:
u = (diff(v)/dt)./v(2:end);
u = u/max(abs(u));

%%
figure;subplot(1,2,1);plot([0,t],v)
xlabel('t (seconds)');
ylabel('I(t)');
title('Input Stimulus')
subplot(1,2,2); plot(t,u)
xlabel('t (seconds)');
ylabel('v(t)');
title('Temporal Contrast')
    
%% Time Encoding
% Set the encoding parameters:
b = 4;          % bias
d = 0.75;       % threshold
C = 0.01;       % capacitance
tauf = 0.015;   % cross-feedback time constant
scale = 0.005;  % amplitude scaling factor

%%
% Encode the signal:
[tk1 tk2]  = iaf_encode_ideal_on_off(u,t,b,d,C,tauf,scale);

%% Time Decoding
% Recover the signal:
u_rec = consistent_decoding_IF_ONOFF(tk1,tk2,t,b,d,C,tauf,scale);

figure;plot(t(tr_vc),u(tr_vc),t(tr_vc),u_rec(tr_vc)); 
xlabel('t (seconds)'); ylabel('u(t)'); 
title('Reconstructed Signal');
legend('original','reconstructed');
  
%%  
% Compute the SNR of the recovered signal:
snr = 10*log10(sum(u(tr_vc).^2)/sum((u(tr_vc)-u_rec(tr_vc)).^2)) 
                                                                 
%%
% _Author: Eftychios A. Pnevmatikakis_
%%
% _Copyright 2009-2010 Trustees of Columbia University_
