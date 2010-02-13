% MIMO_demo
% Performs the example presented in section 4.3 
% The input signal is a vector valued bandlimited stimulus. The stimulus is
% filtered by a kernel that performs arbitrary but known scaling and
% delaying to the inputs. The output of the filter It is encoded with 
% a population of ideal IF neurons. The signal is then reconstructed within 
% the L2 space, without making any bandlimited assumptions

% Author: Eftychios A. Pnevmatikakis
% Copyright 2010 Trustees of Columbia University

% run time ~ 9min on an iMac 5.1, Intel Core 2 Duo, 2.16GHz, 4GB RAM

%% signal generation

dt = 1e-7;
fmax = 100;
t = dt:dt:0.1;
tr_vc = round(0.1*length(t)):round(0.9*length(t)); % truncated vector

mc = floor(0.6*floor(t(end)/dt)*fmax*dt); % # of sinusoidal components

M = 3; % number of inputs
u=zeros(M,length(t));

for i=1:M
    ut = gen_test_signal(t(end)+(2*round(0.1*length(t))-1)*dt,dt,fmax,-Inf,mc);
    u(i,:) = ut((round(0.1*length(t))+1):end-round(0.1*length(t)));
    u(i,:) = u(i,:)/max(abs(u(i,:)));
end

figure;
for i=1:M
    subplot(1,M,i); plot(t,u(i,:));%plot(t_tr,u(i,trun_vect)');
    ylabel('Amplitude');
    xlabel(sprintf('Signal u^{%d}',i));
end


%% neural circuits parameters 

N = 12; % number of neurons
delay = exprnd(1/2/fmax/40,N,M); % delays

scale = 0.5 + 0.5*rand(N,M); % rank should be M

b = 2.3 + rand(1,N);
d = 0.5 + rand(1,N);
C = 0.01*ones(1,N)/2;

%% encoding

% construct dendritic currents
v = zeros(N,length(tr_vc));

for j=1:N
    for i=1:M
        v(j,:) = v(j,:) + scale(j,i)*u(i,tr_vc + round(delay(j,i)/dt));
    end
end    

for j=1:N
    tk = cumsum([0,iaf_encode(v(j,:), dt, b(j), d(j), Inf, C(j))]);
    TK(1:length(tk),j) = tk';
    LN(j) = length(tk);
end

%% decoding
u_rec = consistent_decoding_IF_MIMO(TK, LN, t(tr_vc)-t(tr_vc(1)), b, d, C, N, M, delay, scale);

figure;
for i = 1:M
    subplot(1,M,i);plot(t(tr_vc),u(i,tr_vc),t(tr_vc),u_rec(i,:))
end

for i = 1:M
    snr(i) = 10*log10(sum(u(i,tr_vc).^2)/sum((u(i,tr_vc)-u_rec(i,:)).^2));
end
snr