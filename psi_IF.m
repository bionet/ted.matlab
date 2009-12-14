function psi = psi_IF(tk1,tk2,t)

psi = zeros(1,length(t));

t1 = (tk1 - t).^4;
t2 = (tk2 - t).^4; % (t-tk)^4

fp = find(t>tk2,1);
fm = find(t>tk1,1);

sp=ones(1,length(t));
sm=ones(1,length(t));
sp(fp:end)=-1;
sm(1:fm-1)=-1;
psi = 0.25*(sp.*t2 + sm.*t1);