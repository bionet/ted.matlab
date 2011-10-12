%CONSISTENT_DECODING_IF_ONOFF Decode a signal encoded with an ON-OFF IAF neuron pair.
%   U_REC = CONSISTENT_DECODING_IF_ONOFF(TK1,TK2,T,B,D,C,TAUF,SCALE)
%   decodes a finite energy signal encoded as spike times TK1 and
%   TK2 by a pair of ON-OFF IAF neurons with bias B, threshold D,
%   capacitance C, cross-feedback time constant TAUF, and
%   cross-feedback amplitude scaling factor SCALE.
%
%   The calculation is described in further detail in Section 3.2
%   of the Consistent Recovery paper mentioned in the toolbox
%   references. 

%   Author: Eftychios A. Pnevmatikakis
%   Copyright 2009-2011 Eftychios A. Pnevmatikakis

function u_rec = consistent_decoding_IF_ONOFF(tk1, tk2, t, b, d, C, tauf, scale)

dt = t(2) - t(1);
ln = [length(tk1)-1,length(tk2)-1];
ln2 = cumsum([0,ln]);
q = q_IF_ONOFF(tk1,tk2, b, d, C, tauf, scale, dt);

TK(1:ln(1)+1,1)=tk1';
TK(1:ln(2)+1,2)=tk2';

G = G_pop_IF(TK,ln);
p = p_IF(TK,Inf,ln);
r = r_IF(TK,Inf,ln);

V = [G p r; p' 0 0; r' 0 0];

cv = pinv(V)*[q;0;0];
d0 = cv(end-1);
d1 = cv(end);
c  = cv(1:end-2);

u_rec = d0+d1*t;

for j = 1:2
    for i = 1:ln(j)
        u_rec = u_rec + c(ln2(j)+i)*psi(TK(i,j),TK(i+1,j),t);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% auxiliary functions to generate G, q, p, r, ans psi

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


function q = q_IF_ONOFF(tk1,tk2, b, d, C, tauf, scale, dt)
%Q_IF_ONOFF Compute the t-transform for an ON-OFF IAF neuron pair.
%   Q = Q_IF_ONOFF(TK1,TK2,B,D,C,TAUF,SCALE,DT) computes the
%   t-transform for an ON-OFF IAF neuron pair for the interspike
%   interval [TK1,TK2] using neurons with bias B, threshold D,
%   capacitance C, cross-feedback time constant TAUF, and
%   cross-feedback amplitude scaling factor SCALE.
% 
%   The calculation is described in further detail in Equation 31
%   of the Consistent Recovery paper mentioned in the toolbox
%   references.

q1 = C*d - b*diff(tk1);
fb1 = 0*q1;

for i=1:length(tk1)-1
    ln=length(find(tk2<tk1(i+1)));
    if ln>0
        for j = 1:ln
            fb1(i) = fb1(i) + dt*trapz(cross_fb(max(tk1(i)-tk2(j),0):dt:(tk1(i+1)-tk2(j)),5,7,tauf,scale));
        end
    else
        fb1(i) = 0;
    end
end

q1 = q1 - fb1;

q2 = -C*d + b*diff(tk2);
fb2 = 0*q2;

for i=1:length(tk2)-1
    ln=length(find(tk1<tk2(i+1)));
    if ln>0
        for j = 1:ln
            fb2(i) = fb2(i) + dt*trapz(cross_fb(max(tk2(i)-tk1(j),0):dt:(tk2(i+1)-tk1(j)),5,7,tauf,scale));
        end
    else
        fb2(i) = 0;
    end
end

q2 = q2 + fb2;

q = [q1, q2]';



function h = cross_fb(t,n1,n2,tauf,scale)
%CROSS_FB Compute cross-feedback for ON-OFF IAF neuron pair.
%   H = CROSS_FB(T,N1,N2,TAUF,SCALE) computes the cross-feedback
%   term coupling the two components of an ON-OFF IAF neuron
%   pair over the times T. N1 and N2 denote the exponents for the
%   positive and negative feedback terms, respectively; TAUF
%   denotes the time constant, and SCALE the amplitude of the
%   cross-feedback.
%
%   The calculation is described in further detail in Equation 33 of
%   the Consistent Recovery paper mentioned in the toolbox references.

h = scale/tauf*exp(-t/tauf).*((t/tauf).^n1/factorial(n1)-(t/tauf).^n2/factorial(n2));
