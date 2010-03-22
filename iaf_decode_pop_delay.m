%IAF_DECODE_POP_DELAY Decode a signal encoded by an ensemble of delayed IAF neurons.
%   U_REC = IAF_DECODE_POP_DELAY(TK,LN,T,W,B,D,KD,G,N,DELAY,DT)
%   decodes the signal with bandwidth W rad/s encoded as spike times
%   TK over the times T using N neurons with biases, thresholds,
%   integration constants, and delays respectively specified in the
%   arrays B,D, KD, and DELAY. The number of spikes from each
%   neuron is specified in LN, and the reconstruction matrix is
%   specified in G. The recovered signal is assumed to be
%   sampled at sampling rate 1/DT Hz. 

%   Author: Eftychios A. Pnevmatikakis
%   Copyright 2009-2010 Trustees of Columbia University

function u_rec = iaf_decode_pop_delay(TK,LN,t,W,b,d,kd,G,N,delay,dt)

ln = LN(1:N)-1; % number of interspike intervals of each neuron
ln2 = cumsum([0,ln]);
bwpi = W/pi;

% Perform the t-transform calculation:
for i=1:N 
    q_v(ln2(i)+1:ln2(i+1),1) = ...
    	kd(i)*d(i)-b(i)*(TK(2:LN(i),i)-TK(1:LN(i)-1,i));
end

% Perform the matrix inversion and calculation of coefficients:
ck_v = (G'*pinv(G*G'))*q_v;

% Reconstruct the signal:
u_rec = zeros(1,length(t));

for neur = 1:N
    c = ck_v(ln2(neur)+1:ln2(neur+1),1);
    tsh = (TK(2:LN(neur),neur)+TK(1:LN(neur)-1,neur))/2;
    for i=1:ln(neur)
        u_rec = u_rec + sinc(bwpi*(t-tsh(i)+delay(neur)))*bwpi*c(i);
    end
end
