%% Time Encoding and Decoding with an Integrate-and-Fire Neuron
% This demo illustrates the time encoding and decoding of a
% bandlimited signal using an integrate-and-fire neuron.

%% A Simple Test Signal
% Generate a noiseless signal 0.1 s long sampled at 1 GHz
% containing 3 components no greater than 32 Hz:
dur = 0.1;  
fs = 1e6;
dt = 1/fs;
f = 32; 
bw = 2*pi*f;
t = linspace(0,dur,floor(dur/dt));

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

%% Time Encoding
% The encoding parameters are validated to ensure that signal
% recovery will possible:
b = 3.5;    % bias
d = 0.7;    % threshold
R = 10;     % resistance
C = 0.01;   % capacitance

if ~iaf_recoverable(u,bw,b,d,R,C),
  return
end

fig_title = 'encoding using leaky IAF algorithm';
fprintf(1,'%s\n',fig_title);
s = func_timer(@iaf_encode,u,dt,b,d,R,C);
plot_encoded(t,u,s,fig_title);

%% Time Decoding
% The encoded signal can be recovered using one of several
% different decoding algorithms:
fig_title = 'decoding using leaky IAF algorithm';
fprintf(1,'%s\n',fig_title);
u_rec = func_timer(@iaf_decode,s,dur,dt,bw,b,d,R,C);
plot_compare(t,u,u_rec,fig_title);
%%
M = 5;    % fast decoding parameter

fig_title = 'decoding using leaky fast IAF algorithm';
fprintf(1,'%s\n',fig_title);
u_rec = func_timer(@iaf_decode_fast,s,dur,dt,bw,M,b,d,R,C);
plot_compare(t,u,u_rec,fig_title);
%%

%% Non-leaky Neuron Model
% Setting the neuron's resistance to infinity is equivalent to
% using a non-leaky neuron model for encoding and decoding the
% signal: 
R = inf;

fig_title = 'encoding using nonleaky IAF algorithm';
fprintf(1,'%s\n',fig_title);
s = func_timer(@iaf_encode,u,dt,b,d,R,C);
plot_encoded(t,u,s,fig_title);
%%
fig_title = 'decoding using nonleaky IAF algorithm';
fprintf(1,'%s\n',fig_title);
u_rec = func_timer(@iaf_decode,s,dur,dt,bw,b,d,R,C);
plot_compare(t,u,u_rec,fig_title);
%%
M = 5;    % fast decoding parameter

fig_title = 'decoding using nonleaky fast IAF algorithm';
fprintf(1,'%s\n',fig_title);
u_rec = func_timer(@iaf_decode_fast,s,dur,dt,bw,M,b,d,R,C);
plot_compare(t,u,u_rec,fig_title);

