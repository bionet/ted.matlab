function psi = psi(tk1,tk2,t,varargin)

% psi creates the function reconstruction function psi

% psi(tk1,tk2,t) creates the matrix function reconstruction
% function psi for an ideal IF neuron that corresponds to the to the interspike
% interval [tik1,tk2]. The functions is given by eq. (19).


% psi_LIF(tk1,tk2,t,RC) creates the matrix function reconstruction
% function psi for a LIF neuron with time constant RC that corresponds to 
% the to the interspike interval [tik1,tk2]. The functions is given by eq. (19).

% Inputs
% tk1:  first spike time
% tk2:  second spike time
% t:    time vector
% RC:   time constant of the LIF neuron


% Output
% psi:  the reconstruction function

% Author(s): Eftychios A. Pnevmatikakis
% Copyright 2010 Trustees of Columbia University

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