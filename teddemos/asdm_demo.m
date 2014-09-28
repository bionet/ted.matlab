%% Time Encoding and Decoding with an Asynchronous Sigma-Delta Modulator
% This demo illustrates the time encoding and decoding of a
% bandlimited signal using an asynchronous sigma-delta modulator.

%% Generating a Test Signal
% Generate a noiseless signal 0.1 s long sampled at 1 MHz containing 3
% components no greater than 32 Hz:
dur = 0.1;      % duration
fs = 1e6;       % sampling frequency
dt = 1/fs;      % sampling resolution
f = 32;         
bw = 2*pi*f;    % bandwidth (rad/s)
t = [0:dt:dur]; % time support
np = -inf;      % noise level

if np == -inf,
  fig_title = 'ASDM Input Signal with No Noise';
else
  fig_title = sprintf('ASDM Input Signal with %d dB of Noise',np);
end

rand('twister',0); randn('state',0);
fprintf(1,'%s\n',fig_title);
u = func_timer(@gen_test_signal,dur,dt,f,np);
figure
plot_signal(t,u,fig_title);

%% Time Encoding
% Set the encoding parameters:
b = 3.5;    % bias
d = 0.7;    % threshold
k = 0.01;   % scaling factor

%% 
% Verify that recovery can take place:
if ~asdm_recoverable(u,bw,b,d,k),
  return
end

%%
% Encode the signal:
fig_title = 'Signal Encoded Using ASDM Encoder';
fprintf(1,'%s\n',fig_title);
s = func_timer(@asdm_encode,u,dt,b,d,k);
figure
plot_encoded(t,u,s,fig_title);

%% Time Decoding
% The encoded signal can be recovered using one of several
% different decoding algorithms:
fig_title = 'Signal Decoded Using ASDM Decoder';
fprintf(1,'%s\n',fig_title);
u_rec = func_timer(@asdm_decode,s,dur,dt,bw,b,d,k);
figure
plot_compare(t,u,u_rec,fig_title);
%%
fig_title = 'Signal Decoded Using Threshold-Insensitive ASDM Decoder';
fprintf(1,'%s\n',fig_title);
u_rec = func_timer(@asdm_decode_ins,s,dur,dt,bw,b);
figure
plot_compare(t,u,u_rec,fig_title);
%%
M = 5;    % fast decoding parameter

fig_title = 'Signal Decoded Using Fast ASDM Decoder';
fprintf(1,'%s\n',fig_title);
u_rec = func_timer(@asdm_decode_fast,s,dur,dt,bw,M,b,d,k);
figure
plot_compare(t,u,u_rec,fig_title);

%%
% _Author: Lev Givon_
%%
% _Copyright 2009-2014 Lev Givon_
