% Time encoding and decoding of a bandpass stimulus
% with a gammatone filterbank and an ensemble of 
% integrate-and-fire neurons

% Author: Eftychios A. Pnevmatikakis
% Bionet Group
% Department of Electrical Engineering
% Columbia University
% April 2009

% Note: Malcolm's Slaney auditory toolbox is needed
% http://cobweb.ecn.purdue.edu/~malcolm/interval/1998-010/



clear all;
close all;
addpath AuditoryToolbox

%% sampling rates

dur = 0.25;

Fs = 1000;  % sampling rate of filters
dt_f = 1/Fs;
t_f = dt_f:dt_f:dur;

Ns = 2^9; % oversampling for neural integration
dt = dt_f/Ns; 
Ft = 1/dt;
t = dt:dt:dur;
tr_vc = round(0.1*length(t)):round(0.9*length(t)); % interval of interest

%% Bandpass Signal Construction

fmin = 110;  % minimum frequency of bandpass signal
fmax = 390;  % maximum frequency of bandpass signal   

mc = floor(floor(t(end)/dt)*(fmax-fmin-1)*dt); % maximum sinusoidal components
u =  gen_test_bp_signal(t(end)+(2*round(0.15*length(t))+1)*dt,dt,fmin,fmax,-Inf,mc);
    
u = u((round(0.15*length(t))+1):end-round(0.15*length(t))-1); % truncate first and last values to eliminate discontinuities  

u = u/max(abs(u));  % normalization

figure;subplot(1,2,1);plot(t,u); % plot signal and its fourier transform
        xlabel('Time [sec]'); ylabel('Amplitude'); title('Singal in the time domain');

        
np2 = 2^nextpow2(length(u(tr_vc)));
U = fft(u(tr_vc),np2)/length(u(tr_vc));  % calculate fourier transform

f = Ft/2*linspace(0,1,np2/2+1);
subplot(1,2,2);plot(f,2*abs(U(1:np2/2+1)));
    axis([0,Fs/2,0,1.1*max(2*abs(U(1:np2/2+1)))]);
    xlabel('Frequency [Hz]'); ylabel('Spectrum'); title('Signal in the frequency domain');

%% Filterbank Construction & Filtering

x_n = resample(u,1,Ns);   % resample input signal to the filter rate

Fmin = 100;  % minimum frequency of filterbank
Nf = 16;     % number of filters
Ln = 1.1*length(t_f); 

fcoefs = MakeERBFilters(Fs,Nf,Fmin);  % constract filter coefficients

h = ERBFilterBank([zeros(1,Ln-1) 1 zeros(1,Ln-1)], fcoefs); % calculate impulse responses [eq. (4.4)]
h_resp = 20*log10(abs(fft(h(:,Ln:end)'))); 
freqScale = (0:(Ln-1))/Ln*Fs; 

figure; subplot(1,2,1);plot(freqScale(1:(Ln/2)-1),h_resp(1:Ln/2-1,:)); % plot impulse responses and filterbank gain
    axis([Fmin Fs/2 -60 5]); xlabel('Frequency [Hz]'); ylabel('Filter Response [dB]'); 

if Nf>1
    Gc=sum((10.^(h_resp'/20)).^2);  %filterbank gain
    subplot(1,2,2);plot(freqScale(1:(Ln/2)-1),Gc(1:(Ln/2)-1)); % [Fig. 9]
        xlabel('Frequency [Hz]'); ylabel('Filterbank Gain')
    Gf=mean(Gc(round(1*length(Gc)/8):round(3*length(Gc)/8)));
end

y_n = ERBFilterBank(x_n, fcoefs);  % pass original signal through the filterbank [eq. (2.2)]

%% Synthesis Filters \tilde(h)

t_nr = (-(Ln-1):(Ln-1))*dt_f;
W_l = 2*pi*420;
sinc_n = W_l/pi*sinc(W_l/pi*t_nr);

sinc_f = ERBFilterBank(sinc_n, fcoefs); 

h_bs = sinc_f(:,end:-1:1);  % synthesis filters

%% Time Encoding
y_r = resample(y_n',Ns,1)';

b = logspace(log10(1.3),log10(2.5),Nf); %1.3 + rand(1,Nf);
d = 1 + rand(1,Nf);
kd=0.01*ones(1,Nf);

for i=1:Nf
    tk = [0,cumsum(iaf_encode(y_r(i,:), dt, b(i), kd(i)*d(i)))];%, kd(i));
    %sk = (tk(1:end-1) + tk(2:end))./2; 
    TK(1:length(tk),i)=tk';
    %SK(1:length(sk),i)=sk';
    LN(i)=length(tk);
end

t_f=(-(Ln*Ns-1):(Ln*Ns-1))/Fs/Ns;

ln = LN-1;
ln2= cumsum([0,ln]);

%% t-transform calculation

for i=1:Nf
    tk=TK(1:LN(i),i)';
    q=kd(i)*d(i)-b(i)*diff(tk);
    q_v(ln2(i)+1:ln2(i+1),1)=q;
end

%% G matrix computing  [eq. (3.18)]

H_conv=cell(1,Nf);
for i=1:Nf
    H_conv{1,i} = ERBFilterBank(h_bs(i,:), fcoefs); 
end

HcH = zeros(Nf,Nf,2*Ln-1);

for i=1:Nf
    for j=1:Nf
        hc=H_conv{1,j};
        HcH(i,j,:) = hc(i,:);   %h_i \ast \tilde{h}_j at filter rate
    end
end

for i=1:Nf
    for j=1:Nf
        G(ln2(i)+1:ln2(i+1),ln2(j)+1:ln2(j+1))=G_filterbank(TK(1:LN(i),i)',TK(1:LN(j),j),[zeros(1,Ns-1) squeeze(resample(HcH(i,j,:),Ns,1))],t_f);
    end
end

%% inversion
Ginv = pinv(G,1e-4); %matrix inversion
ck_v = Ginv*q_v;

%% reconstruction

hb_rec=resample(h_bs',Ns,1)';  % reconstruction filter at signal rate

u_temp=zeros(1,size(hb_rec,2));

for i=1:Nf
    tk=TK(1:LN(i),i)';
    sk=(tk(1:end-1)+tk(2:end))/2;
    for j=1:ln(i)
        u_temp = u_temp + ck_v(ln2(i)+j)*circshift(hb_rec(i,:),[0,round(sk(j)*Fs*Ns)]);
                % weight and shift the reconstruction filter at time sk(j)
                % [eq. (3.17)]
    end
end

u_rec = u_temp((Ln-1)*Ns:(Ln-1)*Ns+length(u)-1);

u_rec = u_rec - mean(u_rec(tr_vc));

figure;subplot(1,2,1);plot(t(tr_vc),u(tr_vc),t(tr_vc),u_rec(tr_vc),t(tr_vc),u(tr_vc)-u_rec(tr_vc));
    xlabel('Time [sec]'); ylabel('Amplitude'); title('Recovery in the time domain');
    legend('Original','Recovered','Error');

U_rec = fft(u_rec(tr_vc),np2)/length(u(tr_vc));  % fourier transform of recovered signal

f = Ft/2*linspace(0,1,np2/2+1);
subplot(1,2,2);plot(f,2*abs(U(1:np2/2+1)),f,2*abs(U_rec(1:np2/2+1)),f,abs(2*abs(U(1:np2/2+1))-2*abs(U_rec(1:np2/2+1))));
    axis([0,Fs/2,0,1.1*max(2*abs(U_rec(1:np2/2+1)))]);
    xlabel('Frequency [Hz]'); ylabel('Spectrum'); title('Recovery in the frequency domain');    
    legend('Original','Recovered','Error')
    
%% sequential recovery

ur_i=zeros(Nf,length(u_temp));
ur_s=zeros(Nf,length(u));
for i=1:Nf
    Gi=G(1:ln2(i+1),1:ln2(i+1));
    ck_i = pinv(Gi,1e-4)*q_v(1:ln2(i+1));
    for j=1:i
        tk=TK(1:LN(j),j)';
        sk=(tk(1:end-1)+tk(2:end))/2;
        for k=1:ln(j)
            ur_i(i,:)=ur_i(i,:)+ck_i(ln2(j)+k)*circshift(hb_rec(j,:),[0,round(sk(k)*Fs*Ns)]);
        end
    end
    ur_s(i,:)=ur_i(i,(Ln-1)*Ns:(Ln-1)*Ns+length(u)-1);
    ur_s(i,:) = ur_s(i,:) - mean(ur_s(i,tr_vc));
end


%% Plot all 16 recovered signals

figure;
for i=1:Nf
    subplot(4,Nf/4,i);plot(t(tr_vc),u(tr_vc),t(tr_vc),ur_s(i,tr_vc))
        xlabel('Time [sec]'); ylabel('Amplitude');
        title(sprintf('# of Neurons: %d',i));
        if i==1
            legend('Original','Recovered');
        end
end

%% Plot recovered signals using 1,2,3,4,8 and 16 neurons
figure;
    subplot(2,3,1);plot(t(tr_vc),u(tr_vc),t(tr_vc),ur_s(1,tr_vc));
        %xlabel('Time [sec]'); 
        title('# of Neurons: 1'); ylabel('Amplitude');
    subplot(2,3,2);plot(t(tr_vc),u(tr_vc),t(tr_vc),ur_s(2,tr_vc));
        %xlabel('Time [sec]'); 
        title('# of Neurons: 2');
    subplot(2,3,3);plot(t(tr_vc),u(tr_vc),t(tr_vc),ur_s(3,tr_vc));
        %xlabel('Time [sec]'); 
        title('# of Neurons: 3');
    subplot(2,3,4);plot(t(tr_vc),u(tr_vc),t(tr_vc),ur_s(4,tr_vc));
        xlabel('Time [sec]'); title('# of Neurons: 4'); ylabel('Amplitude');
    subplot(2,3,5);plot(t(tr_vc),u(tr_vc),t(tr_vc),ur_s(8,tr_vc));
        xlabel('Time [sec]'); title('# of Neurons: 8');
    subplot(2,3,6);plot(t(tr_vc),u(tr_vc),t(tr_vc),ur_s(16,tr_vc));
        xlabel('Time [sec]'); title('# of Neurons: 16');

%% MSE and SNR

tr_vt = round(0.1*length(u)):round(0.9*length(u));

for i=1:Nf
    ms=u(tr_vt)-ur_s(i,tr_vt);
    mse(i)=10*log10(mean(ms.^2));
    snr(i)=10*log10(sum(u(tr_vt).^2)/sum(ms.^2));
end

figure;plot(1:Nf,mse); grid on;
xlabel('# of Neurons'); ylabel('MSE (dB)');
title('MSE as a Function of the Number of Neurons');

figure;plot(1:Nf,snr); grid on;
xlabel('# of Neurons'); ylabel('SNR (dB)');
title('SNR as a Function of the Number of Neurons');