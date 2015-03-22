%GEN_DIRICHLET_COEFFS Generate random Dirichlet coefficients for a real signal
%   AM = GEN_DIRICHLET_COEFFS(M) generates 2*M+1 random Dirichlet
%   coefficients for a real signal. The first M coefficients are
%   the complex conjugates of the last M coefficients reversed.

%   Author: Lev Givon
%   Copyright 2009-2015 Lev Givon

function am = gen_dirichlet_coeffs(M)

am = rand(1,2*M+1)+j*rand(1,2*M+1); 
am(1:M) = conj(am(end:-1:end-M+1));
am(M+1) = rand(1);
