%% Time Encoding and Decoding with Multiple IAF Neurons
% This demo illustrates the time encoding and decoding of a
% bandlimited signal using an ensemble of integrate-and-fire 
% neurons.

%% A Simple Test Signal
% Generate a noiseless signal 0.1 s long sampled at 1 GHz
% containing 3 components no greater than 32 Hz:
dur = 0.1;  
fs = 1e6;
dt = 1/fs;
f = 32; 
bw = 2*pi*f;
t = [0:dt:dur];

np = -inf;   

if np == -inf,
  fig_title = 'IAF input signal with no noise';
else
  fig_title = sprintf('IAF input signal with %d dB of noise',np);
end
rand('twister',0); randn('state',0);
fprintf(1,'%s\n',fig_title);
u = func_timer(@gen_test_signal,dur,dt,f,np);
plot_signal(t,u,fig_title);

%% Encoding and Decoding with Leaky Neurons
% In this example, the input signal is encoded using two leaky IAF
% neurons with different encoding parameters. The parameters are
% validated to ensure that signal recovery will be possible:
b1 = 3.5;    % bias
d1 = 0.7;    % threshold
R1 = 10;     % resistance
C1 = 0.01;   % capacitance

if ~iaf_recoverable(u,bw,b1,d1,R1,C1),
  return
end

fig_title = 'encoding using leaky IAF algorithm (encoder #1)';
fprintf(1,'%s\n',fig_title);
s1 = func_timer(@iaf_encode,u,dt,b1,d1,R1,C1);
plot_encoded(t,u,s1,fig_title);
%%
b2 = 3.4;    % bias
d2 = 0.8;    % threshold
R2 = 9.0;    % resistance
C2 = 0.01;   % capacitance

if ~iaf_recoverable(u,bw,b2,d2,R2,C2),
  return
end

fig_title = 'encoding using leaky IAF algorithm (encoder #2)';
fprintf(1,'%s\n',fig_title);
s2 = func_timer(@iaf_encode,u,dt,b2,d2,R2,C2);
plot_encoded(t,u,s2,fig_title);
%%
fig_title = 'decoding using leaky IAF algorithm';
fprintf(1,'%s\n',fig_title);
u_rec = func_timer(@iaf_decode_pop,{s1,s2},dur,dt,bw, ...
                  {b1,b2},{d1,d2},{R1,R2},{C1,C2});
plot_compare(t,u,u_rec,fig_title);

%% Encoding and Decoding with Non-leaky Neurons
% The above example is repeated with neurons whose resistance is
% infinite:
R1 = inf;
R2 = inf;
%%
fig_title = 'encoding using non-leaky IAF algorithm (encoder #1)';
fprintf(1,'%s\n',fig_title);
s1 = func_timer(@iaf_encode,u,dt,b1,d1,R1,C1);
plot_encoded(t,u,s1,fig_title);
%%
fig_title = 'encoding using non-leaky IAF algorithm (encoder #2)';
fprintf(1,'%s\n',fig_title);
s2 = func_timer(@iaf_encode,u,dt,b2,d2,R2,C2);
plot_encoded(t,u,s2,fig_title);
%%
fig_title = 'decoding using non-leaky IAF algorithm';
fprintf(1,'%s\n',fig_title);
u_rec = func_timer(@iaf_decode_pop,{s1,s2},dur,dt,bw, ...
                   {b1,b2},{d1,d2},{R1,R2},{C1,C2});
plot_compare(t,u,u_rec,fig_title);
%%

