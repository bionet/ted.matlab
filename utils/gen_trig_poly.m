%GEN_TRIG_POLY Generate a trigonometric polynomial test signal
%   U = GEN_TRIG_POLY(T,DT,M) generates a trigonometric
%   polynomial of length T s with 2*M+1 frequency
%   components sampled at time resolution DT s.
%
%   U = GEN_TRIG_POLY(T,DT,AM) generates a trigonometric polynomial
%   using the 2*M+1 Dirichlet coefficients in AM.

%   Author: Lev Givon
%   Copyright 2009-2012 Lev Givon

function u = gen_trig_poly(T,dt,am)

if isscalar(am),
    M = am;
    am = gen_dirichlet_coeffs(M);
else
    N = length(am);
    if mod(N,2) == 0,
        error('number of coefficients must be odd');
    end
    M = floor(N/2);
end

if M < 1,
    error('number of coefficients must be at least 1');
end

u_fft = zeros(1,floor(T/dt));
u_fft(end-M+1:end) = am(1:M)*sqrt(T)/dt;
u_fft(1:M+1) = am(M+1:2*M+1)*sqrt(T)/dt;

u = real(ifft(u_fft));
