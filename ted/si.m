%SI Sine integral function.
%   Y = SI(X) computes the sine integral INT(SIN(T)/T,T,0,X) using
%   a series approximation. This implementation is faster than 
%   Matlab's SININT function.
%
%   See also SININT.

%   Adapted from specfun (http://www.netlib.org/specfun/)
%   by Lev Givon
%   Copyright 2009-2011 Lev Givon

function y = si(x)

y = arrayfun(@si_scalar, x);

function y = si_scalar(x);

sgn = sign(x);
x = abs(x);

p2 = 1.570796326794897;
el = 0.5772156649015329;
eps = 1e-15;

bj = zeros(101);
x2 = x*x;

if x == 0.0,
    ci = inf;
    si = 0.0;
elseif x <= 16.0,
    xr = -0.25*x2;
    ci = el+log(x)+xr;
    for k=2:40,
        xr = -0.5*xr*(k-1)/(k*k*(2*k-1))*x2;
        ci = ci+xr;
        if abs(xr) < abs(ci)*eps,
            break
        end
    end
    xr = x;
    si = x;
    for k=1:40,
        xr = -0.5*xr*(2*k-1)/k/(4*k*k+4*k+1)*x2;
        si = si+xr;
        if abs(xr) < abs(si)*eps,
            y = sgn*si;
            return
        end
    end
elseif x <= 32.0,
    m = fix(47.2+0.82*x);
    xa1 = 0.0;
    xa0 = 1.0e-100;
    for k=m:-1:1,
        xa = 4.0*k*xa0/x-xa1;
        bj(k) = xa;
        xa1 = xa0;
        xa0 = xa;
    end
    xs = bj(1);
    for k=3:2:m,
        xs = xs+2.0*bj(k);
    end
    bj(1) = bj(1)/xs;
    for k=2:m,
        bj(k) = bj(k)/xs;
    end
    xr = 1.0;
    xg1 = bj(1);
    for k=2:m,
        xr = 0.25*xr*(2.0*k-3.0)^2/((k-1.0)*(2.0*k-1.0)^2)*x;
        xg1 = xg1+bj(k)*xr;
    end
    xr = 1.0;
    xg2 = bj(1);
    for k=2:m,
        xr = 0.25*xr*(2.0*k-5.0)^2/((k-1.0)*(2.0*k-3.0)^2)*x;
        xg2 = xg2+bj(k)*xr;
    end
    xcs = cos(x/2.0);
    xss = sin(x/2.0);
    ci = el+log(x)-x*xss*xg1+2*xcs*xg2-2*xcs*xcs;
    si = x*xcs*xg1+2*xss*xg2-sin(x);
else
    xr = 1.0;
    xf = 1.0;
    for k=1:9,
        xr = -2.0*xr*k*(2*k-1)/x2;
        xf = xf+xr;
    end
    xr = 1.0/x;
    xg = xr;
    for k=1:8,
        xr = -2.0*xr*(2*k+1)*k/x2;
        xg = xg+xr;
    end
    ci = xf*sin(x)/x-xg*cos(x)/x;
    si = p2-xf*cos(x)/x-xg*sin(x)/x;
end
y = si*sgn;
