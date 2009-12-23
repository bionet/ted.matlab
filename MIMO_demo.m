%% Shifting + Scaling Example MIMO Consistent Recovery
clear all
close all

%% signal generation

dt = 1e-7;
fmax = 100;
Omega = 2*pi*fmax;
T = pi/Omega;
t = dt:dt:20*T;
sc = 0.6;

M = 3; % number of inputs
u=zeros(M,length(t));

for i=1:M
    u(i,:) = gen_trig_signal(20*T,dt,fmax,sc);
end
%clear A; clear B;

trun_vect=round((dt:dt:18*T)/dt);  %truncated vector
t_tr=t(trun_vect); %truncated time

figure;
for i=1:M
    subplot(1,M,i); plot(t,u(i,:));%plot(t_tr,u(i,trun_vect)');
    ylabel('Amplitude');
    xlabel(sprintf('Signal u^{%d}',i));
    %axis([min(t_tr) max(t_tr) 0.5*floor(min(min(u))/0.5)
    %0.5*ceil(max(max(u))/0.5)])
end


%% neural circuits parameters 

N = 9; % number of neurons
delay = exprnd(T/20,N,M); % delays

scale = 0.5 + 0.5*rand(N,M); % rank should be M

b = 2.3 + rand(1,N);
d = 0.5 + rand(1,N);
C = 0.01*ones(1,N)/2;

%% encoding

% construct dendritic currents
v = zeros(N,length(trun_vect));

for j=1:N
    for i=1:M
        v(j,:) = v(j,:) + scale(j,i)*u(i,trun_vect + round(delay(j,i)/dt));
    end
end    

for j=1:N
    tk = iaf_encode_perfect(v(j,:), t(trun_vect)-t(trun_vect(1)), b(j), d(j), C(j));
    TK(1:length(tk),j) = tk';
    LN(j) = length(tk);
end

%% decoding
u_rec = consistent_decoding_IF_MIMO(TK, LN, t_tr-t(trun_vect(1)), b, d, C, N, M, delay, scale);

figure;
for i = 1:M
    subplot(1,M,i);plot(t,u(i,:),t_tr,u_rec(i,:))
end