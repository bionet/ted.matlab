function u_rec = LIF_pop_decode_S2(s_list, dur, dt, b, delta, R, C, n, lamda, N)
% LIF_pop_decode_S1 reconstruct signals encoded with an ensemble of LIF
% neuron with random thresholds
%
% u_rec = LIF_pop_decode_S1(TK, LN, b, d, R, C, lamda, N) reconstructs the encoded 
% stimulus u that belongs in the Sobolev space S2 and is encoded with a
% population on LIF neurons with random thresholds. The process is described in detail
% in section 4.2.2
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
% Copyright 2009-2012 Eftychios A. Pnevmatikakis

t = [0:dt:dur];
RC = R.*C;


ts_list = cellfun(@cumsum,s_list,'UniformOutput',false);
LN = cellfun(@length,ts_list);
TK = zeros(max(LN),N);

for i = 1:N
    TK(1:LN(i),i) = ts_list{i}.';
end


if sum(n.^2) == 0 % deal with the case when noise is absent
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

G = G_S2_pop_LIF(TK,LN,R,C,n,N);
F = F_S2_pop_LIF(TK,LN,R,C,n,N);

%% determine coefficients

c = zeros(sum(ln),Nl); 
d = zeros(2,Nl);

[Q R] = qr(F);
Q1 = Q(:,1:2);
Q2 = Q(:,3:end);
R  = R(1:2,:);

for l = 1:Nl
    M = G + sum(ln)*lamda(l)*eye(sum(ln));
    S = Q2'*M*Q2;
    Sinv = S'*pinv(S*S');
    c(:,l) = Q2*Sinv*Q2'*q_v;
    d(:,l) = inv(R)*Q1'*(q_v-M*c(:,l));    
end

%% reconstruction

u_rec = diag(d(1,:))*ones(Nl,length(t)) + d(2,:)'*t;

for i = 1:N
    for j=1:ln(i)
        psi = psi_S2_LIF(TK(j,i),TK(j+1,i),RC(i),t,dt)/nC(i);
        u_rec = u_rec + c(ln2(i)+j,:)'*psi;
    end
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% auxiliary functions to generate G, F, q and psi


function G = G_S2_pop_LIF(TK,LN,R,C,n,N)
% G_S2_pop_LIF create the matrix G [equation (47)]
%
% G = G_S2_pop_LIF(tk,RC) creates the matrix G with entries G^{ij}[k,l] =
% <phi^i_k,psi^j_l> for the reconstruction of a stimulus that belongs in the
% Sobolev space S2 and is encoded with a population of leaky integrate-and-fire neurons.
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
%
% Author(s): Eftychios A. Pnevmatikakis, based on a SciPy implementation by
%               Robert J. Turetsky


ln = LN(1:N)-1;
ln2 = [0 cumsum(ln)];
RC = R.*C;
nC = n.*C;

G = zeros(sum(ln),sum(ln));
G_diag = G; % diagonal entries of the block matrix G^{ii}

for i = 1:N
    G_diag(ln2(i)+1:ln2(i+1),ln2(i)+1:ln2(i+1)) = G_S2_LIF(TK(1:LN(i),i)',RC(i))/(nC(i)^2);
end

GK  = @(ta,tb,RCi,RCj,tip,tjp)(RCi*RCj*exp(-(tip-ta)/RCi - (tjp-tb)/RCj));
GFi = @(ta,tb,RCi,RCj,tip,tjp)(((RCi-RCj+tb)*(RCi*ta-RCi^2-ta^2/2)+ta^3/6)*GK(ta,tb,RCi,RCj,tip,tjp));
GFj = @(ta,tb,RCi,RCj,tip,tjp)(((RCj-RCi+ta)*(RCj*tb-RCj^2-tb^2/2)+tb^3/6)*GK(ta,tb,RCi,RCj,tip,tjp));


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
                    Gij(k,l) = -GFj(tim,tjm,RCi,RCj,tip,tjp) + GFj(tim,tjp,RCi,RCj,tip,tjp) + GFj(tip,tjm,RCi,RCj,tip,tjp) ...
                                -GFj(tip,tjp,RCi,RCj,tip,tjp);
                elseif (tjm<=tim)+(tim<=tjp)+(tjp<=tip)==3
                    Gij(k,l) = -GFj(tim,tjm,RCi,RCj,tip,tjp) + GFi(tim,tjp,RCi,RCj,tip,tjp) + GFj(tip,tjm,RCi,RCj,tip,tjp) ...
                                -GFj(tip,tjp,RCi,RCj,tip,tjp) + GK(tjp,tjp,RCi,RCj,tip,tjp)/(RCi+RCj)*(RCi^4) ...
                                -GK(tim,tim,RCi,RCj,tip,tjp)/(RCi+RCj)*(RCj^4);        
                elseif (tjm<=tim)+(tip<tjp)==2
                    Gij(k,l) = -GFj(tim,tjm,RCi,RCj,tip,tjp) + GFi(tim,tjp,RCi,RCj,tip,tjp) + GFj(tip,tjm,RCi,RCj,tip,tjp) ...
                                -GFi(tip,tjp,RCi,RCj,tip,tjp) + GK(tip,tip,RCi,RCj,tip,tjp)/(RCi+RCj)*(RCj^4) ...
                                -GK(tim,tim,RCi,RCj,tip,tjp)/(RCi+RCj)*(RCj^4);  
                elseif (tim<=tjm)+(tjp<=tip)==2
                    Gij(k,l) = -GFi(tim,tjm,RCi,RCj,tip,tjp) + GFi(tim,tjp,RCi,RCj,tip,tjp) + GFj(tip,tjm,RCi,RCj,tip,tjp) ...
                                -GFj(tip,tjp,RCi,RCj,tip,tjp) + GK(tjp,tjp,RCi,RCj,tip,tjp)/(RCi+RCj)*(RCi^4) ...
                                -GK(tjm,tjm,RCi,RCj,tip,tjp)/(RCi+RCj)*(RCi^4);                                   
                elseif (tim<=tjm)+(tjm<=tip)+(tip<=tjp)==3
                    Gij(k,l) = -GFi(tim,tjm,RCi,RCj,tip,tjp) + GFi(tim,tjp,RCi,RCj,tip,tjp) + GFj(tip,tjm,RCi,RCj,tip,tjp) ...
                                -GFi(tip,tjp,RCi,RCj,tip,tjp) + GK(tip,tip,RCi,RCj,tip,tjp)/(RCi+RCj)*(RCj^4) ...
                                -GK(tjm,tjm,RCi,RCj,tip,tjp)/(RCi+RCj)*(RCi^4);                
                else
                    Gij(k,l) = -GFi(tim,tjm,RCi,RCj,tip,tjp) + GFi(tim,tjp,RCi,RCj,tip,tjp) + GFi(tip,tjm,RCi,RCj,tip,tjp) ...
                                -GFi(tip,tjp,RCi,RCj,tip,tjp);        
                end
                %Gij(k,l) = Gij(k,l)*exp(-tip/RCi-tjp/RCj);
            end
        end
        G(ln2(i)+1:ln2(i+1),ln2(j)+1:ln2(j+1)) = Gij/nC(i)/nC(j);
    end
end

G = G + G' + G_diag;


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


function F = F_S2_pop_LIF(TK,LN,R,C,n,N)
% F_S1_LIF create the matrix F [equation (47)]
%
% F = F_S1_LIF(tk,RC) creates the vector F with entries F^i[k,1] =
% <phi^i_k,1> and F^i[k,2] = <phi^i_k,t> for the reconstruction of a stimulus that belongs in the
% Sobolev space S2 and is encoded with a population of leaky integrate-and-fire neurons.
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
F = zeros(sum(ln),2);

for i = 1:N
    F(ln2(i)+1:ln2(i+1),:) = F_S2_LIF(TK(1:LN(i),i)',RC(i))/nC(i);
end

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



