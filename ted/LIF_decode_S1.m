function u_rec = LIF_decode_S1(s, dur, dt, b, delta, R, C, lamda)
% LIF_decode_S1 reconstruct signals encoded with a LIF neuron with random threshold
%
% u_rec = LIF_decode_S1(t, b, d, R, C, lamda) reconstructs the encoded 
% stimulus u that belongs in the Sobolev space S1 and is encoded with a 
% LIF neuron with random threshold. The process is described in detail
% in section 3.2.1
%
% Inputs:
% s:      inter-spike intervals
% dur:    duration of recovery
% dt:     reverse of sample rate
% delta:  mean threshold of the LIF neuron
% b:      bias
% R:      resistance 
% C:      capacitance
% lamda:  smoothing parameter [can be a vector]
%
% Output:
% u_rec:  reconstructed stimulus. Each row corresponds to the reconstructed
%         function for a different choice of the smoothing parameter  
%
% Author: Eftychios A. Pnevmatikakis
% Copyright 2009-2011 Eftychios A. Pnevmatikakis

t = [0:dt:dur];
RC = R*C;
Nl = length(lamda);
tk = cumsum(s);
ln = length(tk)-1; % # of samples

q = q_LIF(tk,b, delta, R, C); % t-transform

G = G_S1_LIF(tk,RC);    % G-matrix
F = F_S1_LIF(tk,RC);    % F-matrix

c = zeros(ln,Nl); 
d = zeros(1,Nl);

for l = 1:Nl
    M = G+ln*lamda(l)*eye(ln);
    Minv = M'*pinv(M*M');
    S = F'*Minv*F;
    Sinv = S'*pinv(S*S');
    d(l) = Sinv*F'*Minv*q;
    c(:,l) = Minv*(eye(ln)-F*Sinv*F'*Minv)*q;
end

u_rec = diag(d)*ones(Nl,length(t));

for i=1:ln
    psi = psi_S1_LIF(tk(i),tk(i+1),RC,t,dt);
    u_rec = u_rec + c(i,:)'*psi;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% auxiliary functions to generate G, F, q and psi

function G = G_S1_LIF(tk,RC)
% G_S1_LIF create the matrix G [equation (28)]
%
% G = G_S1_LIF(tk,RC) creates the matrix G with entries G[i,j] =
% <phi_i,psi_j> for the reconstruction of a stimulus that belongs in the
% Sobolev space S1 and is encoded with a leaky integrate-and-fire neuron.
%
% Inputs
% tk:  vector of spike times
% RC:  time constant of the LIF neuron (R*C)
%
% Output
% G: the matrix G


ln = length(tk)-1;

G = zeros(ln,ln);

for i=1:ln-1
    for j=i+1:ln
        G(i,j)=RC*(1-exp(-(tk(j+1)-tk(j))/RC))*(RC*(tk(i+1)-tk(i)*exp(-(tk(i+1)-tk(i))/RC)) - ...
            RC^2*(1-exp(-(tk(i+1)-tk(i))/RC)));
    end
end

G = G + G';

for i=1:ln
    G(i,i) = RC^2*tk(i)*(1-exp(-(tk(i+1)-tk(i))/RC))^2 + RC^2*(tk(i+1)-tk(i))...
        - 2*RC^3*(1-exp(-(tk(i+1)-tk(i))/RC)) + 0.5*RC^3*(1-exp(-2*(tk(i+1)-tk(i))/RC));
end


function F = F_S1_LIF(tk,RC)
% F_S1_LIF create the matrix F [equation (28)]
%
% F = F_S1_LIF(tk,RC) creates the vector F with entries F[i] =
% <phi_i,1> for the reconstruction of a stimulus that belongs in the
% Sobolev space S1 and is encoded with a leaky integrate-and-fire neuron.
% Inputs
%
% tk:  vector of spike times
% RC:  time constant of the LIF neuron (R*C)
%
% Output
% F: the vector F


F=RC*(1 - exp(-diff(tk)/RC))';



function q = q_LIF(tk,b,d,R,C)
% q_LIF create the t-transform for a LIF neuron [equation (8)]
%
% q = q_LIF(tk,b,d,R,C) creates the vector q with entries the right
% hand side of the t-transform of a signal encoded with a LIF neuron.
%
% Inputs
% tk:  vector of spike times
% b :  bias
% d :  threshold
% R :  resistance
% C :  capacitacne
%
% Output
% q :  the vector q

q = (C*(d+b*R*(exp(-diff(tk)/(R*C))-1)))';



function psi = psi_S1_LIF(tk1,tk2,RC,t,dt)
% ps1_S1_LIF create the stimulus representation function [equation (27)]
%
% psi = psi_S1_LIF(tk1,tk2,RC,t,dt) creates the stimulus representation
% function for the reconstruction of a stimulus that belongs in the
% Sobolev space S1 and is encoded with a leaky integrate-and-fire neuron.
% The representation function corresponds to sample that is obtained from 
% the interspike interval [tk1,tk2]
%
% Inputs
% tk1: time of first spike
% tk2: time of second spike
% RC:  time constant of the LIF neuron (R*C)
% t:   time vector [0,1]
% dt:  time step
%
% Output
% psi: the representation function

% split time interval [0,1] into [0,tk1],(tk1,tk2] and (tk2,1]
ind1 = zeros(1,length(t));
ind2 = zeros(1,length(t));
ind3 = zeros(1,length(t));

ind1(1:round(tk1/dt)) = 1;
ind2(round(tk1/dt)+1:round(tk2/dt)) = 1;
ind3(round(tk2/dt)+1:end) = 1;

% assign corresponding values
a1 = ind1.*t*RC*(1 - exp(-(tk2-tk1)/RC));
a2 = ind2.*(RC*t - RC^2*exp(-(tk2-t)/RC) + RC*(RC-tk1)*exp(-(tk2-tk1)/RC));
a3 = ind3.*(RC*(tk2-tk1*exp(-(tk2-tk1)/RC)) - RC^2*(1-exp(-(tk2-tk1)/RC)));

% concatenate
psi = a1 + a2 + a3;
