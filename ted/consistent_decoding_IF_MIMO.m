%CONSISTENT_DECODING_IF_MIMO Decode several signal encoded with an ensemble of IAF neurons.
%   U_REC = CONSISTENT_DECODING_IF_MIMO(TK,LN,T,B,D,C,N,M,DELAY,SCALE)
%   decodes a vector-valued signal comprising M inputs encoded as
%   spike times TK over the times T using N ideal IAF neurons with
%   biases, thresholds, and capacitances respectively specified in the
%   arrays B, D, and C and a filtering kernel that delays and scales
%   the inputs using the values in the DELAY and SCALE matrices.
%   delaying and scaling to the inputs. The number of spikes from each
%   neuron is specified in LN. 
%
%   The calculation is described in further detail in Section 4.2
%   of the Consistent Recovery paper mentioned in the toolbox
%   references.

%   Author: Eftychios A. Pnevmatikakis
%   Copyright 2009-2011 Eftychios A. Pnevmatikakis

function u_rec = consistent_decoding_IF_MIMO(TK,LN,t,b,d,C,N,M,delay,scale)

dt = t(2) - t(1);
ln = LN-1;
ln2 = cumsum([0,ln]);

q_v = zeros(ln2(end),1);

for i=1:N
    q_v(ln2(i)+1:ln2(i+1),1) = C(i)*d(i) - b(i)*diff(TK(1:LN(i),i)');    
end

G = zeros(ln2(end),ln2(end));
p = zeros(ln2(end),M);
r = zeros(ln2(end),M);

for inp = 1:M
    G = G + G_pop_IF(TK-repmat(delay(:,inp)',max(LN),1),ln,scale(:,inp));
    p(:,inp) = p_IF(TK-repmat(delay(:,inp)',max(LN),1),Inf,ln,scale(:,inp))'; 
    r(:,inp) = r_IF(TK-repmat(delay(:,inp)',max(LN),1),Inf,ln,scale(:,inp))';
end

V = [G p r; p' zeros(M,2*M); r' zeros(M,2*M)];

cv = pinv(V)*[q_v;zeros(2*M,1)];
d0 = cv(ln2(end)+1:ln2(end)+M);
d1 = cv(ln2(end)+M+1:end);
c  = cv(1:ln2(end));

u_rec = zeros(M,length(t));
u_rec = repmat(d0,1,length(t)) + d1*t;

for j = 1:N
    for spk = 1:ln(j)
        for inp = 1:M
            u_rec(inp,:) = u_rec(inp,:) + scale(j,inp)*c(ln2(j)+spk)*psi(TK(spk,j)-delay(j,inp),TK(spk+1,j)-delay(j,inp),t);
        end
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% auxiliary functions to generate G, p, r, ans psi

function G = G_pop_IF(TK,ln,varargin)
%G_POP_IF Compute the reconstruction matrix G for multiple IAF neuron decoder.
%   G = G_POP_IF(TK,LN,W) computes the matrix G with entries 
%   G[i,j] = <phi_k^i,psi_l^j> used to reconstruct a signal in L2
%   space that was encoded with a population of ideal IAF
%   neurons. TK contains the spike times, while LN denotes the
%   number of spikes due to each neuron. If specified, W denotes
%   the weights on the inputs prior to encoding; if not specified,
%   the weights are all assumed to be 1.
%
%   The calculation is described in further detail in Equation 29 
%   of the Consistent Recovery paper mentioned in the toolbox 
%   references.

N = length(ln);
G = zeros(sum(ln),sum(ln));
ln2 = cumsum([0,ln]);

if nargin > 2
    w = varargin{1};
else
    w = ones(1,N);
end

for i = 1:N
    for j = 1:N
        Gb = G_block_IF(TK(1:ln(i)+1,i)',TK(1:ln(j)+1,j)');
        G(ln2(i)+1:ln2(i+1),ln2(j)+1:ln2(j+1)) = w(i)*w(j)*Gb;
    end
end



function Gb = G_block_IF(ti,tj)
%G_BLOCK_IF Compute the reconstruction matrix for multiple IAF time decoder.
%   G = G_BLOCK_IF(TI,TJ) computes the reconstruction matrix
%   G[i,j] = <phi_k^i,psi_l^j> used to decode a signal in L2 space
%   that was encoded by a population of ideal IAF neurons. TI and
%   TJ contain the spike times from neurons i and j, respectively.
%
%   The calculation is described in further detail in Equation 29 of the
%   Consistent Recovery paper mentioned in the toolbox references.

li=length(ti)-1;
lj=length(tj)-1;

Gb = zeros(li,lj);

if isequal(ti,tj)
    for i=1:li
        for j=1:li
            tmz=ti(min(i,j));
            tmp=ti(min(i,j)+1);
            tpz=ti(max(i,j));
            tpp=ti(max(i,j)+1);
            Gb(i,j)=((tpp-tmz)^5-(tpp-tmp)^5+(tpz-tmp)^5-(tpz-tmz)^5)/20;
        end
    end
else
    for k=1:li
        for l=1:lj
            tip=ti(k+1);
            tim=ti(k);
            tjp=tj(l+1);
            tjm=tj(l);
            if tjp<=tim
                Gb(k,l)=((tip-tjm)^5+(tim-tjp)^5-(tim-tjm)^5-(tip-tjp)^5)/20;
            elseif (tjm<=tim)+(tim<=tjp)+(tjp<=tip)==3
                Gb(k,l)=((tip-tjm)^5-(tip-tjp)^5+(tjp-tim)^5-(tim-tjm)^5)/20;
            elseif (tjm<=tim)+(tip<tjp)==2
                Gb(k,l)=((tip-tjp)^5+(tip-tjm)^5-(tim-tjp)^5-(tim-tjm)^5)/20;
            elseif (tim<=tjm)+(tjm<=tip)==2
                Gb(k,l)=((tip-tjm)^5+(tip-tjm)^5-(tip-tjp)^5-(tim-tjp)^5)/20;
            elseif (tim<=tjm)+(tjm<=tip)+(tip<=tjp)==3
                Gb(k,l)=((tip-tjp)^5+(tip-tjm)^5-(tim-tjp)^5+(tim-tjm)^5)/20;
            else
                Gb(k,l)=((tip-tjp)^5+(tim-tjm)^5-(tim-tjp)^5-(tip-tjm)^5)/20;
            end
        end
    end
end



function p = p_IF(TK,RC,varargin)
%P_IF Compute the inner product <1,\phi_k>.
%   P = P_IF(TK,RC) computes the inner product <1,\phi_k> for
%   a leaky IAF neuron with time constant RC (an ideal neuron may
%   be simulated by setting RC to Inf). The spike times generated
%   by the neuron are specified in TK.
%
%   P = P_IF(TK,RC,LN,W) computes the inner product for a population
%   of IAF neurons; LN contains the number of spikes from each
%   neuron, and W contains the values by which each neuron weights
%   the input. If W is not specified, the input is not weighted.


if nargin > 2
    ln = varargin{1};
else
    ln = length(TK)-1;
    w = 1;
    N = 1;
end
if nargin > 3
    w = varargin{2};
    N = length(ln);
else
    w = ones(1,length(ln));
    N = length(ln);
end
if nargin > 4
    error('Too many input arguments.');
end

ln2 = cumsum([0,ln]);
p = zeros(sum(ln),1);

if isinf(RC)
    for i = 1:N
        p(ln2(i)+1:ln2(i+1)) = w(i)*diff(TK(1:ln(i)+1,i));
    end
else
    for i = 1:N
        p(ln2(i)+1:ln2(i+1)) = w(i)*RC(i)*(1-exp(-diff(TK(1:ln(i)+1,i))/RC(i)));
    end
end




function r = r_IF(TK,RC,varargin)

%R_IF Compute the inner product <t,\phi_k>.
%   R = R_IF(TK,RC) computes the inner product <t,\phi_k> for a
%   leaky IAF neuron with time constant RC (an ideal neuron may be
%   simulated by setting RC to Inf). The spike times generated by
%   the neuron are specified in TK.
%
%   R = R_IF(TK,RC,LN,W) computes the inner product for a
%   population of IAF neurons; LN contains the number of spikes
%   from each neuron, and W contains the values by which each
%   neuron weights the input. If W is not specified, the input is
%   not weighted.

if nargin > 2
    ln = varargin{1};
else
    ln = length(TK)-1;
    w = 1;
    N = 1;
end
if nargin > 3
    w = varargin{2};
    N = length(ln);
else
    w = ones(1,length(ln));
    N = length(ln);
end
if nargin > 4
    error('Too many input arguments.');
end

ln2 = cumsum([0,ln]);
r = zeros(sum(ln),1);

if isinf(RC)
    for i = 1:N
        r(ln2(i)+1:ln2(i+1)) = w(i)*diff(TK(1:ln(i)+1,i).^2)/2;
    end
else
    for i = 1:N
        r(ln2(i)+1:ln2(i+1)) = w(i)*(RC(i)^2)*(TK(2:ln(i)+1,i)/RC(i)-1-(TK(1:ln(i),i)/RC(i)-1).*exp(-diff(TK(1:ln(i)+1,i))/RC(i)));
    end
end



function psi = psi(tk1,tk2,t,varargin)

%PSI Create the signal reconstruction function psi.
%   P = PSI(TK1,TK2,T,RC) computes the matrix reconstruction function
%   psi for an IAF neuron for the interspike interval [TK1,TK2] over
%   the times T. If the neuron's time constant RC is not specified, it
%   is assumed to be infinite.
%   
%   The calculation is described in further detail in Equation 19 of
%   the Consistent Recovery paper mentioned in the toolbox references.

psi = zeros(1,length(t));

if nargin > 3
    RC = varargin{1};
elseif nargin == 3
    RC = Inf;
else
    error('Wrong number of input arguments.');
end


if isinf(RC)
    t1 = (tk1 - t).^4;
    t2 = (tk2 - t).^4; 

    fp = find(t>tk2,1);
    fm = find(t>tk1,1);

    sp = ones(1,length(t));
    sm = ones(1,length(t));
    sp(fp:end)=-1;
    sm(1:fm-1)=-1;
    psi = 0.25*(sp.*t2 + sm.*t1);

else    
    
    t1 = (tk1 - t)/RC;
    t2 = (tk2 - t)/RC;

    f=@(x)(x.^3-3*x.^2+6*x-6);

    ft  = (RC^4)*f([t1;t2]);
    etk = exp(-(tk2-tk1)/RC);
    ett = 12*(RC^4)*exp(-[t1;t2]);

    fp = find(t>tk2,1);
    fm = find(t>tk1,1);

    psi(1:fm-1)  = ft(2,1:fm-1)-ft(1,1:fm-1)*etk;
    psi(fm:fp-1) = ett(2,fm:fp-1)+ft(2,fm:fp-1)+ft(1,fm:fp-1)*etk;
    psi(fp:end)  = ft(1,fp:end)*etk-ft(2,fp:end);
end
