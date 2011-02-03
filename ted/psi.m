%PSI Create the signal reconstruction function psi.
%   P = PSI(TK1,TK2,T,RC) computes the matrix reconstruction function
%   psi for an IAF neuron for the interspike interval [TK1,TK2] over
%   the times T. If the neuron's time constant RC is not specified, it
%   is assumed to be infinite.
%   
%   The calculation is described in further detail in Equation 19 of
%   the Consistent Recovery paper mentioned in the toolbox references.

%   Author: Eftychios A. Pnevmatikakis
%   Copyright 2009-2011 Eftychios A. Pnevmatikakis

function psi = psi(tk1,tk2,t,varargin)

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
