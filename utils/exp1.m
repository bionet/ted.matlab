%EXP1 Exponential integral function.
%   Y = EXP1(Z) computes the exponential integral E1 
%   of a complex value.

%   Adapted from specfun (http://www.netlib.org/specfun/ 
%   by Lev Givon
%   Copyright 2009-2014 Lev Givon

function y = exp1(z)

y = arrayfun(@exp1_scalar, z);

function y = exp1_scalar(z)

el = 0.5772156649015328;
x = real(z);
a0 = abs(z);

if a0 == 0,
    ce1 = complex(inf, 0);
elseif (a0 <= 10.0) || ((x < 0.0) && (a0 < 20.0)),
    ce1 = complex(1.0, 0.0);
    cr = complex(1.0, 0.0);
    for k=1:150,
        cr = -cr*k*z/(k+1.0)^2;
        ce1 = ce1+cr;
        if abs(cr) <= abs(ce1)*1e-15,
            break;
        end
    end
    ce1 = -el-log(z)+z*ce1;
else
    ct0 = complex(0.0, 0.0);
    for k=120:-1:1,
        ct0 = k/(1.0+k/(z+ct0));
    end
    ct = 1.0/(z+ct0);
    ce1 = exp(-z)*ct;
    if (x <= 0.0) && (imag(z) == 0.0),
        ce1 = ce1-pi*complex(0.0, 1.0);
    end
end
if isnan(ce1),
    y = nan;
else
    y = ce1;
end

