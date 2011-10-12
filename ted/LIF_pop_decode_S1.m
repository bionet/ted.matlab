function u_rec = LIF_pop_decode_S1(s_list, dur, dt, b, delta, R, C, n, lamda, N)

% LIF_pop_decode_S1 reconstruct signals encoded with an ensemble of LIF
% neuron with random thresholds
%
% u_rec = LIF_pop_decode_S1(TK, LN, b, d, R, C, lamda, N) reconstructs the encoded 
% stimulus u that belongs in the Sobolev space S1 and is encoded with a
% population on LIF neurons with random thresholds. The process is described in detail
% in section 4.2.1
%
% Inputs:
% s_list  cell array of inter-spike intervals, one cell for one neuron
% dur:    duration of recovery
% dt:     reverse of sample rate
% delta:  mean threshold of the LIF neuron
% b:      bias
% R:      resistance 
% C:      capacitance
% n:      noise power
% lamda:  smoothing parameter [can be a vector]
% N:      number of neurons
%
% Output:
% u_rec:  reconstructed stimulus. Each row corresponds to the reconstructed
%         function for a different choice of the smoothing parameter  
%
% Author: Eftychios A. Pnevmatikakis
% Copyright 2009-2011 Eftychios A. Pnevmatikakis

t = [0:dt:dur];
RC = R.*C;


ts_list = cellfun(@cumsum,s_list,'UniformOutput',false);
LN = cellfun(@length,ts_list);
TK = zeros(max(LN),N);

for i = 1:N
    TK(1:LN(i),i) = ts_list{i}.';
end


if sum(n.^2) == 0
    n = ones(1,N);
end
nC = n.*C;

ln = LN(1:N)-1;
ln2 = [0 cumsum(ln)];

Nl = length(lamda);


%% t-transform
for i = 1:N
    q = q_LIF(TK(1:LN(i),i)',b(i),delta(i),R(i),C(i))/nC(i);
    q_v(ln2(i)+1:ln2(i+1),1) = q;    
end

%% create G and F
G = G_S1_pop_LIF(TK,LN,R,C,n,N);
F = F_S1_pop_LIF(TK,LN,R,C,n,N);

%% determine coefficients

c = zeros(sum(ln),Nl); 
d = zeros(1,Nl);

for l = 1:Nl
    M = G + sum(ln)*lamda(l)*eye(sum(ln));
    Minv = M'*pinv(M*M');
    S = F'*Minv*F;
    Sinv = S'*pinv(S*S');
    d(l) = Sinv*F'*Minv*q_v;
    c(:,l) = Minv*(eye(sum(ln))-F*Sinv*F'*Minv)*q_v;
end

%% reconstruction

u_rec = diag(d)*ones(Nl,length(t));

for i = 1:N
    for j=1:ln(i)
        psi = psi_S1_LIF(TK(j,i),TK(j+1,i),RC(i),t,dt)/nC(i);
        u_rec = u_rec + c(ln2(i)+j,:)'*psi;
    end
end






%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% auxiliary functions to generate G, F, q and psi

function G = G_S1_pop_LIF(TK,LN,R,C,n,N)

% G_S1_pop_LIF create the matrix G [equation (47)]
%
% G = G_S1_pop_LIF(tk,RC) creates the matrix G with entries G^{ij}[k,l] =
% <phi^i_k,psi^j_l> for the reconstruction of a stimulus that belongs in the
% Sobolev space S1 and is encoded with a population of leaky integrate-and-fire neurons.
%
% Inputs
% TK:  vector of spike trains
% LN:  number of spikes per neuron
% R:   resistance
% C:   capacitance
% n:   noise power per neuron
% N:   number of neurons
%
% Output
% G: the matrix G


ln = LN(1:N)-1;
ln2 = [0 cumsum(ln)];
RC = R.*C;
nC = n.*C;

G = zeros(sum(ln),sum(ln));
G_diag = G; % diagonal entries of the block matrix G^{ii}

for i = 1:N
    G_diag(ln2(i)+1:ln2(i+1),ln2(i)+1:ln2(i+1)) = G_S1_LIF(TK(1:LN(i),i)',RC(i))/(nC(i)^2);
end

for i = 2:N
    ti = TK(1:LN(i),i)';
    RCi = RC(i);
    li = ln(i);
    for j = 1:i-1;
        Gij = zeros(ln(i),ln(j));
        tj = TK(1:LN(j),j)';
        RCj = RC(j);
        lj = ln(j);
        RCij = RCi*RCj/(RCi+RCj);
        for k=1:li
            for l=1:lj
                tip=ti(k+1);
                tim=ti(k);
                tjp=tj(l+1);
                tjm=tj(l);
                if tjp<=tim
                    Gij(k,l) = RCi*RCj*(exp(tip/RCi)-exp(tim/RCi))*(exp(tjp/RCj)*(tjp-RCj)-exp(tjm/RCj)*(tjm-RCj));
                elseif (tjm<=tim)+(tim<=tjp)+(tjp<=tip)==3
                    Gij(k,l) = RCi*RCj*(exp(tip/RCi)-exp(tim/RCi))*(exp(tim/RCj)*(tim-RCj)-exp(tjm/RCj)*(tjm-RCj)) ...
                           + RCj*(RCij*(exp(tjp/RCij)*(tjp-RCij)-exp(tim/RCij)*(tim-RCij))-RCj*RCij*(exp(tjp/RCij)-exp(tim/RCij)) - ...
                                    RCi*exp(tim/RCj)*(tim-RCj)*(exp(tjp/RCi)-exp(tim/RCi))) ...
                           + RCi*(RCij*(exp(tjp/RCij)*(tjp-RCij)-exp(tim/RCij)*(tim-RCij))-RCi*RCij*(exp(tjp/RCij)-exp(tim/RCij)) - ...
                                    RCj*exp(tim/RCi)*(tim-RCi)*(exp(tjp/RCj)-exp(tim/RCj))) ...
                           + RCi*RCj*(exp(tip/RCi)-exp(tjp/RCi))*(exp(tjp/RCj)*(tjp-RCj)-exp(tim/RCj)*(tim-RCj));         
                elseif (tjm<=tim)+(tip<tjp)==2
                    Gij(k,l) = RCi*RCj*(exp(tip/RCi)-exp(tim/RCi))*(exp(tim/RCj)*(tim-RCj)-exp(tjm/RCj)*(tjm-RCj)) ...
                           + RCj*(RCij*(exp(tip/RCij)*(tip-RCij)-exp(tim/RCij)*(tim-RCij))-RCj*RCij*(exp(tip/RCij)-exp(tim/RCij)) - ...
                                    RCi*exp(tim/RCj)*(tim-RCj)*(exp(tip/RCi)-exp(tim/RCi))) ...
                           + RCi*(RCij*(exp(tip/RCij)*(tip-RCij)-exp(tim/RCij)*(tim-RCij))-RCi*RCij*(exp(tip/RCij)-exp(tim/RCij)) - ...
                                    RCj*exp(tim/RCi)*(tim-RCi)*(exp(tip/RCj)-exp(tim/RCj))) ...
                           + RCi*RCj*(exp(tjp/RCj)-exp(tip/RCj))*(exp(tip/RCi)*(tip-RCi)-exp(tim/RCi)*(tim-RCi)); 
                elseif (tim<=tjm)+(tjp<=tip)==2
                    Gij(k,l) = RCi*RCj*(exp(tjp/RCj)-exp(tjm/RCj))*(exp(tjm/RCi)*(tjm-RCi)-exp(tim/RCi)*(tim-RCi)) ...
                           + RCi*(RCij*(exp(tjp/RCij)*(tjp-RCij)-exp(tjm/RCij)*(tjm-RCij))-RCi*RCij*(exp(tjp/RCij)-exp(tjm/RCij)) - ...
                                    RCj*exp(tjm/RCi)*(tjm-RCi)*(exp(tjp/RCj)-exp(tjm/RCj))) ...
                           + RCj*(RCij*(exp(tjp/RCij)*(tjp-RCij)-exp(tjm/RCij)*(tjm-RCij))-RCj*RCij*(exp(tjp/RCij)-exp(tjm/RCij)) - ...
                                    RCi*exp(tjm/RCj)*(tjm-RCj)*(exp(tjp/RCi)-exp(tjm/RCi))) ...
                           + RCi*RCj*(exp(tip/RCi)-exp(tjp/RCi))*(exp(tjp/RCj)*(tjp-RCj)-exp(tjm/RCj)*(tjm-RCj));                 
                elseif (tim<=tjm)+(tjm<=tip)+(tip<=tjp)==3
                    Gij(k,l) = RCi*RCj*(exp(tjp/RCj)-exp(tjm/RCj))*(exp(tjm/RCi)*(tjm-RCi)-exp(tim/RCi)*(tim-RCi)) ...
                           + RCi*(RCij*(exp(tip/RCij)*(tip-RCij)-exp(tjm/RCij)*(tjm-RCij))-RCi*RCij*(exp(tip/RCij)-exp(tjm/RCij)) - ...
                                    RCj*exp(tjm/RCi)*(tjm-RCi)*(exp(tip/RCj)-exp(tjm/RCj))) ...
                           + RCj*(RCij*(exp(tip/RCij)*(tip-RCij)-exp(tjm/RCij)*(tjm-RCij))-RCj*RCij*(exp(tip/RCij)-exp(tjm/RCij)) - ...
                                    RCi*exp(tjm/RCj)*(tjm-RCj)*(exp(tip/RCi)-exp(tjm/RCi))) ...
                           + RCi*RCj*(exp(tjp/RCj)-exp(tip/RCj))*(exp(tip/RCi)*(tip-RCi)-exp(tjm/RCi)*(tjm-RCi));                
                else
                    Gij(k,l) = RCi*RCj*(exp(tjp/RCj)-exp(tjm/RCj))*(exp(tip/RCi)*(tip-RCi)-exp(tim/RCi)*(tim-RCi));                
                end
                Gij(k,l) = Gij(k,l)*exp(-tip/RCi-tjp/RCj);
            end
        end
        G(ln2(i)+1:ln2(i+1),ln2(j)+1:ln2(j+1)) = Gij/nC(i)/nC(j);
    end
end

G = G + G' + G_diag;


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




function F = F_S1_pop_LIF(TK,LN,R,C,n,N)
% F_S1_LIF create the matrix F [equation (47)]
%
% F = F_S1_LIF(tk,RC) creates the vector F with entries F^i[k] =
% <phi^i_k,1> for the reconstruction of a stimulus that belongs in the
% Sobolev space S1 and is encoded with a population of leaky integrate-and-fire neurons.
%
% Inputs
% TK:  vector of spike trains
% LN:  number of spikes per neuron
% R:   resistance
% C:   capacitance
% n:   noise power per neuron
% N:   number of neurons
%
% Output
% F: the vector F


RC = R.*C;
nC = n.*C;
ln = LN(1:N) - 1;
ln2 = [0 cumsum(ln)];
F = zeros(sum(ln),1);

for i = 1:N
    F(ln2(i)+1:ln2(i+1),1) = F_S1_LIF(TK(1:LN(i),i)',RC(i))/nC(i);
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
