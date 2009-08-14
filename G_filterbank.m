function G=G_filterbank2(tk1,tk2,HH,t)

% Construction of the G^{ij} block of the G matrix in the
%    example of encoding with a gammatone filterbank:
% Inputs
% tk1: Spike train from neuron i
% tk2: Spike train from neuron j
% HH: vector of h^i \ast \tilde{h}^j
% t: appropriate time vector of HH

% Output
% G: G^{ij} block

% Author: Eftychios A. Pnevmatikakis
% Bionet Group
% Department of Electrical Engineering
% Columbia University
% April 2009

dt=t(2)-t(1);
sk1 = (tk1(1:end-1) + tk1(2:end))/2;
sk2 = (tk2(1:end-1) + tk2(2:end))/2;

Hc=dt*cumtrapz(HH);

G = zeros(length(sk1), length(sk2));
for j=1:length(sk1)
    for k=1:length(sk2)
        a1=find(t<(tk1(j)-sk2(k)),1,'last');
        a2=find(t<(tk1(j+1)-sk2(k)),1,'last');
        G(j,k) = Hc(a2)-Hc(a1-1);
    end
end