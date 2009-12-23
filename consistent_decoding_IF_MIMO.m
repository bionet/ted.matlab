function u_rec = consistent_decoding_IF_MIMO(TK, LN, t, b, d, C, N, M, delay, scale)

% TK spike times of each neuron
% LN number of spikes for each neuron
% t time vector
% b bias vector
% d thresholds vector 
% C capacitance vector
% N number of neurons
% M number of inputs
% delay matrix of delays
% scale matrix of scalings
% dt time step

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
    G = G + G_pop_IF_scale(TK+repmat(delay(:,inp)',max(LN),1),ln,scale(:,inp),N);
    p(:,inp) = p_pop_IF_scale(TK+repmat(delay(:,inp)',max(LN),1),ln,scale(:,inp),N);
    r(:,inp) = r_pop_IF_scale(TK+repmat(delay(:,inp)',max(LN),1),ln,scale(:,inp),N);
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
            u_rec(inp,:) = u_rec(inp,:) + scale(j,inp)*c(ln2(j)+spk)*psi_IF(TK(spk,j)+delay(j,inp),TK(spk+1,j)+delay(j,inp),t);
        end
    end
end