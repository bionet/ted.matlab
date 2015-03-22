%% Time Encoding and Decoding Using Gammatone Filters and the Trigonometric Polynomial Approximation

%% 
% Generate a trigonometric polynomial input signal:
M = 200;
Omega = 2*pi*1200;
T = 2*pi*M/Omega;
dur = T;
dt = 1e-7;
t = [0:dt:dur-dt];
rand('twister',0); randn('state',0); 
u_orig = gen_trig_poly(T,dt,M);
u = u_orig/max(abs(u_orig));

%% 
% Generate and display the filter bank:
N = 16;
h = make_gammatone_fb(t,N); % size(h) == [N, length(t)]

%%
% Filter the input signal with the filter bank:
v = zeros(size(h));
for n=1:N,
    v(n,:) = filter_trig_poly(u,h(n,:));    
end

%% 
% Encode the filtered signal:
b = 2+0.1*rand(1,N);
d = 0.001*ones(1,N);
s = {};
for n=1:N,
    s{end+1} = iaf_encode(v(n,:),dt,b(n),d(n));
end

%% 
% Decode the filtered signal with different numbers of neurons.
% To save time decoding, the encoded signal is recovered at a
% coarse time resolution and then upsampled to the original resolution:
n_neurons = [2, 6, 10, 16];
dt_temp = dt*100;
t_temp = [0:dt_temp:dur-dt_temp];
for n=n_neurons,
    u_rec_temp = iaf_decode_filt_trig_pop(s(1:n),dur,dt_temp,Omega,M,b(1:n),d(1:n),h);
    u_rec = interp1(t_temp,u_rec_temp,t,'spline');
    figure();
    plot_compare(t,u,u_rec,sprintf('Recovery with %i Neurons',n));
end

%%
% _Author: Lev Givon_
%%
% _Copyright 2009-2015 Lev Givon_
