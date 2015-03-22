%IAF_DECODE_FILT_TRIG_POP Decode a signal encoded with an ensemble of IAF neurons.  
%   U_REC = IAF_DECODE_FILT_TRIG_POP(S,DUR,DT,OMEGA,M,B,D,H,R,C)
%   decodes the signal with duration DUR s and bandwidth OMEGA rad/s
%   encoded as a cell array of spike interval sequences S using
%   integrate-and-fire neurons with biases B, firing thresholds D,
%   resistances R, and capacitances C. The recovered signal is assumed
%   to have been sampled at sampling rate 1/DT Hz and convolved with
%   the impulse functions in the rows of H. This function assumes that
%   the signal can be approximated as a sum of 2*M+1 trigonometric
%   polynomials.
%
%   The parameter H is optional; if specified, the encoded signal is
%   assumed to have been filtered by the impulse responses in the
%   rows of H before being encoded.
%
%   The parameters R and C are also optional. If all of the
%   resistances are inf (the default), a non-leaky neuron model is
%   used. The capacitances are all assumed to be 1 if not specified.

%   Author: Lev Givon
%   Copyright 2009-2015 Lev Givon

function u_rec = iaf_decode_filt_trig_pop(s,dur,dt,Omega,M,b,d,varargin)

% Number of interspike interval trains:
N = length(s);

R = inf(1,N);
C = ones(1,N);
if nargin >= 8,
    h = varargin{1};
end
if nargin >= 9,
    R = varargin{2};
end
if nargin >= 10,
    C = varargin{3};
end
if nargin >= 11,
    error('Too many input arguments.');
end

TM = 2*pi*M/Omega;
OmegaM = Omega/M;
sqrtTM = sqrt(TM);
em = @(m,t) exp(1j*m*OmegaM*t)/sqrtTM;

F = [];
q = [];
if all(isinf(R)), % ideal case
    for i=1:N,

        % Get the Dirichlet coefficients of the filter:
        if exist('h'),
            hm = get_dirichlet_coeffs(h(i,:),dt,Omega,M);
        end
        ts = cumsum(s{i});
        K = length(s{i})-1;
        F_temp = zeros(K,2*M+1);
        for k=1:K;
            for m=-M:M,
                if m == 0,
                    F_temp(k,m+M+1) = s{i}(k+1)/sqrtTM;
                else
                    F_temp(k,m+M+1) = (em(m,ts(k+1))-em(m,ts(k)))/(1j*m*OmegaM);
                end
            end                
            if exist('hm'),
                F_temp(k,:) = (hm*sqrtTM/dt).*F_temp(k,:);
            end
        end
        F = [F; F_temp];
        q = [q, C(i)*d(i)-b(i)*s{i}(2:end)];
    end
else            % leaky case
    for i=1:N,

        % Get the Dirichlet coefficients of the filter:
        if exist('h'),
            hm = get_dirichlet_coeffs(h(i,:),dt,Omega,M);
        end
        ts = cumsum(s{i});
        K = length(s{i})-1;
        F_temp = zeros(K,2*M+1);
        RC = R(i)*C(i);
        for k=1:K;
            for m=-M:M,
                if m == 0,
                    F_temp(k,m+M+1) = (exp(ts(k+1)/RC)-exp(ts(k)/RC))*...
                        exp(-ts(k+1)/RC)/(sqrtTM/RC);
                else
                    x = 1j*m*OmegaM+1/RC;
                    F_temp(k,m+M+1) = (exp(ts(k+1)*x)-exp(ts(k)*x))*...
                        exp(-ts(k+1)/RC)/(sqrtTM*x);
                end
            end                
            if exist('hm'),
                F_temp(k,:) = (hm*sqrtTM/dt).*F_temp(k,:);
            end
        end
        F = [F; F_temp];
        q = [q, C(i)*d(i)-b(i)*RC*(1-exp(-s{i}(2:end)/RC))];
    end
end

c = pinv(F)*q';
t = [0:dt:dur-dt];
u_rec = zeros(1,length(t));
for m=-M:M,
    u_rec = u_rec + c(m+M+1)*em(m,t);
end
u_rec = real(u_rec);
