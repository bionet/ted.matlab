function u_rec = LIF_decode_S2(s, dur, dt, b, delta, R, C, lamda)
% (tk,t,dt,lamda,R,C,delta,b)
% LIF_decode_S2 reconstruct signals encoded with a LIF neuron with random threshold
%
% u_rec = LIF_decode_S2(t, b, d, R, C, lamda) reconstructs the encoded 
% stimulus u that belongs in the Sobolev space S2 and is encoded with a 
% LIF neuron with random threshold. The process is described in detail
% in section 3.2.2
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

G = G_S2_LIF(tk,RC);    % G-matrix
F = F_S2_LIF(tk,RC);    % F-matrix

[Q R] = qr(F);      % QR Decomposition
Q1 = Q(:,1:2);
Q2 = Q(:,3:end);
R  = R(1:2,:);

c = zeros(ln,Nl); 
d = zeros(2,Nl);

for l = 1:Nl
    M = G+ln*lamda(l)*eye(ln);
    S = Q2'*M*Q2;
    Sinv = S'*pinv(S*S');
    c(:,l) = Q2*Sinv*Q2'*q;
    d(:,l) = inv(R)*Q1'*(q-M*c(:,l));
end

u_rec = diag(d(1,:))*ones(Nl,length(t)) + d(2,:)'*t;

for i=1:ln
    psi = psi_S2_LIF(tk(i),tk(i+1),RC,t,dt);
    u_rec = u_rec + c(i,:)'*psi;
end    





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% auxiliary functions to generate G, F, q and psi

function G = G_S2_LIF(tk,RC)
% G_S2_LIF create the matrix G [equation (39)]
%
% G = G_S2_LIF(tk,RC) creates the matrix G with entries G[i,j] =
% <phi_i,psi_j> for the reconstruction of a stimulus that belongs in the
% Sobolev space S2 and is encoded with a leaky integrate-and-fire neuron.
%
% Inputs
% tk:  vector of spike times
% RC:  time constant of the LIF neuron (R*C)
%
% Output
% G: the matrix G

ln = length(tk)-1;

yk = 1 - exp(-diff(tk)/RC);
zk = tk(2:end) - RC - (tk(1:end-1)-RC).*exp(-diff(tk)/RC);
wk = diff(tk)/RC;


G = zeros(ln,ln);

for i=1:ln-1
    for j=i+1:ln
        G(i,j) = 1/3*tk(i)^3*yk(i)*yk(j) - 1/2*tk(i)^2*(yk(i)*zk(j) + yk(j)*zk(i))+ ...
            tk(i)*zk(i)*zk(j)+...
        zk(j)*(((tk(i+1)-tk(i))/2-RC)*(tk(i+1)-tk(i))+RC^2*yk(i))-...
        yk(j)*((tk(i+1)-RC)*(tk(i+1)^2-tk(i)^2)/2-(tk(i+1)^3-tk(i)^3)/3 + RC^2*zk(i));
    end
end

G = G + G';

for i=1:ln
     G(i,i) =1/3*(tk(i)^3)*(yk(i)^2) - (tk(i)^2)*yk(i)*zk(i) + tk(i)*(zk(i))^2 + ((tk(i+1)-tk(i))^3)/3 ...
         - RC*(tk(i+1)-tk(i))^2 - (RC^2)*(1-2*yk(i))*(tk(i+1)-tk(i))+...
         0.5*(RC^3)*(1-exp(-2*(tk(i+1)-tk(i))/RC));
end

G = RC^2*G;



function F = F_S2_LIF(tk,RC)
% F_S2_LIF create the matrix F [equation (36)]
%
% F = F_S1_LIF(tk,RC) creates the matrix F with entries F[i,1] =
% <phi_i,1> and F[i,2] = <phi_i,t> for the reconstruction of 
% a stimulus that belongs in the Sobolev space S1 and is encoded with a 
% leaky integrate-and-fire neuron.
%
% Inputs
% tk:  vector of spike times
% RC:  time constant of the LIF neuron (R*C)
%
% Output
% F: the matrix F

ln = length(tk)-1;

F = zeros(ln,2);
F(:,1)=RC*(1 - exp(-diff(tk)/RC))';
F(:,2)=RC*(tk(2:end)-RC - (tk(1:end-1)-RC).*exp(-(diff(tk))/RC))';



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


function psi = psi_S2_LIF(tk1,tk2,RC,t,dt)
% ps1_S2_LIF create the stimulus representation function [equation (35)]
%
% psi = psi_S2_LIF(tk1,tk2,RC,t,dt) creates the stimulus representation
% function for the reconstruction of a stimulus that belongs in the
% Sobolev space S2 and is encoded with a leaky integrate-and-fire neuron.
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
a1 = ind1.*RC.*(-1/6*t.^3*(1 - exp(-(tk2-tk1)/RC))+1/2*t.^2*(tk2-RC - (tk1-RC)*exp(-(tk2-tk1)/RC)));

a2 = RC*ind2.*(t/2.*((RC^2+(t-RC).^2).*exp(-(tk2-t)/RC)-(RC^2+(tk1-RC)^2)*exp(-(tk2-tk1)/RC)) - ...
    1/6*(((t-RC).^3+RC^2*(3*t-5*RC)).*exp(-(tk2-t)/RC)-((tk1-RC)^3+RC^2*(3*tk1-5*RC))*exp(-(tk2-tk1)/RC))-...
    1/6*(t.^3).*(1 - exp(-(tk2-t)/RC))+1/2*(t.^2).*(tk2-RC - (t-RC).*exp(-(tk2-t)/RC)));

a3 = ind3.*RC.*(t/2*(RC^2+(tk2-RC)^2-(RC^2+(tk1-RC)^2)*exp(-(tk2-tk1)/RC)) - ...
    1/6*((tk2-RC)^3+RC^2*(3*tk2-5*RC)-((tk1-RC)^3+RC^2*(3*tk1-5*RC))*exp(-(tk2-tk1)/RC)));

% concatenate
psi = a1 + a2 + a3;
