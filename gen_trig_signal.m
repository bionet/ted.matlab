function u = gen_trig_signal(dur,dt,fmax,sc)

% generate a trigonometric function
% dur = signal duration
% dt = time resolution
% fmax = bandwidth [Hz]
% sc = sparsity coefficient [0-1]

M = ceil(dur*fmax);
dur_new = M/fmax;

Nc = ceil(M*sc);

aj = rand(1,Nc)-0.5;
ph = 2*pi*rand(1,Nc);
ic = randperm(M);

t = dt:dt:dur_new;
u = 0*t;
wm = 2*pi*fmax/M;

for i = 1:Nc
    u = u + 2*aj(i)*cos(ic(i)*wm*t+ph(i));
end

u = u(1:round(dur/dt));
u = u/max(abs(u));