%% Time Encoding and Decoding with an Integrate-and-Fire Neuron
% This demo illustrates the time encoding and decoding of a
% bandlimited signal using an integrate-and-fire neuron.

%% Generating a Test Signal
% Generate a noiseless signal 0.1 s long sampled at 1 GHz containing 3
% components no greater than 32 Hz:
dur = 0.1;      % duration
fs = 1e6;       % sampling frequency
dt = 1/fs;      % sampling resolution
f = 32;         
bw = 2*pi*f;    % bandwidth (rad/s)
t = [0:dt:dur]; % time support
np = -inf;      % noise level

if np == -inf,
  fig_title = 'IAF input signal with no noise';
else
  fig_title = sprintf('IAF input signal with %d dB of noise',np);
end

rand('twister',0); randn('state',0);
fprintf(1,'%s\n',fig_title);
u = func_timer(@gen_test_signal,dur,dt,f,np);
plot_signal(t,u,fig_title);

%% Time Encoding
% The IAF time encoder can make use of a leaky or non-leaky neuron
% model (i.e., when the neuron's resistance is infinite). Both models
% are demonstrated below.
%
% Set the encoding parameters:
b = 3.5;  % bias
d = 0.7;  % threshold
R = 10;   % resistance
C = 0.01; % capacitance

%%
% Verify that recovery can take place with the leaky and non-leaky
% parameters:
if ~iaf_recoverable(u,bw,b,d,R,C),
  return
end

if ~iaf_recoverable(u,bw,b,d,inf,C),
  return
end

%%
% Encode the signal using the leaky model: 
fig_title = 'encoding using leaky IAF algorithm';
fprintf(1,'%s\n',fig_title);
s_leaky = func_timer(@iaf_encode,u,dt,b,d,R,C);
plot_encoded(t,u,s_leaky,fig_title);

%%
% Encode the signal using the non-leaky model:
fig_title = 'encoding using nonleaky IAF algorithm';
fprintf(1,'%s\n',fig_title);
s_nonleaky = func_timer(@iaf_encode,u,dt,b,d,inf,C);
plot_encoded(t,u,s_nonleaky,fig_title);

%% Time Decoding
% The signal can be recovered for both the leaky and non-leaky models:
fig_title = 'decoding using leaky IAF algorithm';
fprintf(1,'%s\n',fig_title);
u_rec_leaky = func_timer(@iaf_decode,s_leaky,dur,dt,bw,b,d,R,C);
plot_compare(t,u,u_rec_leaky,fig_title);
%%
fig_title = 'decoding using nonleaky IAF algorithm';
fprintf(1,'%s\n',fig_title);
u_rec_nonleaky = func_timer(@iaf_decode,s_nonleaky,dur,dt,bw,b,d,inf,C);
plot_compare(t,u,u_rec_nonleaky,fig_title);

%%
% Decoding can also be performed using a faster algorithm:
M = 5; % fast decoding parameter

fig_title = 'decoding using leaky fast IAF algorithm';
fprintf(1,'%s\n',fig_title);
u_rec_leaky = func_timer(@iaf_decode_fast,s_leaky,dur,dt,bw,M,b,d,inf,C);
plot_compare(t,u,u_rec_leaky,fig_title);
%%
fig_title = 'decoding using nonleaky fast IAF algorithm';
fprintf(1,'%s\n',fig_title);
u_rec_nonleaky = func_timer(@iaf_decode_fast,s_nonleaky,dur,dt,bw,M,b,d,inf,C);
plot_compare(t,u,u_rec_nonleaky,fig_title);

%%
% _Author: Lev Givon_
%%
% _Copyright 2009-2010 Trustees of Columbia University_
