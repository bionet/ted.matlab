%GET_DIRICHLET_COEFFS Compute Dirichlet coefficients of a signal
%   AM = GET_DIRICHLET_COEFFS(U,DT,OMEGA,M) computes the M Dirichlet
%   coefficients AM of a trigonometric polynomial with bandwidth OMEGA
%   rad/s using the FFT given one entire period U of the
%   polynomial. The time resolution of the signal is assumed to be DT
%   s. 

%   Author: Lev Givon
%   Copyright 2009-2012 Lev Givon

function am = get_dirichlet_coeffs(u,dt,Omega,M)

T = 2*pi*M/Omega;
u_fft = fft(u);
am = zeros(1,2*M+1);
am(1:M) = u_fft(end-M+1:end)*dt/sqrt(T);
am(M+1:2*M+1) = u_fft(1:M+1)*dt/sqrt(T);
