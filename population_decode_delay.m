% Recover a signal encoded by a population of ideal IAF neurons with
% delay.
%
% Inputs
%   TK    - spike times of all neurons
%   LN    - number of spikes of each neuron
%   t     - time vector
%   W     - bandwidth
%   b     - biases
%   d     - thresholds
%   kd    - integration constants
%   N     - number of neurons used for recovery
%   delay - delays of each neuron
%   dt    - time step
%
% Output
%   u_rec - recovered signal

% Author: Eftychios A. Pnevmatikakis
% Copyright 2009 Trustees of Columbia University

function u_rec = population_decode_delay(TK, LN, t, W, b, d, kd, G, ...
                                         N, delay, dt)

ln = LN(1:N)-1; % number of interspike intervals of each neuron
ln2 = cumsum([0,ln]);
bwpi = W/pi;

% Perform the t-transform calculation:
for i=1:N 
    q_v(ln2(i)+1:ln2(i+1),1) = kd(i)*d(i)-b(i)*(TK(2:LN(i),i)-TK(1:LN(i)-1,i));
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