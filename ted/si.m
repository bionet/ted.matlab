%SI Sine integral function.
%  Y = SI(X) computes the sine integral INT(SIN(T)/T,T,0,X).
%  This implementation is faster than Matlab's SININT function
%
%  See also SININT.

%  Adapted from specfun (http://www.netlib.org/specfun/)
%  by Lev Givon
%  Copyright 2009 Trustees of Columbia University

function y = si(x)

y = arrayfun(@si_scalar, x);

function y = si_scalar(x);

x2 = x*x;
if x == 0,
    y = 0.0;
else
    if x <= 1.0,
        y = (((((3.1e-7)*x2-2.834e-5)*x2+1.66667e-3)*x2-5.555556e-2)*x2+1.0)*x;
    else
        fx = ((((x2+38.027264)*x2+265.187033)*x2+335.67732)*x2+38.102495)/ ...
             ((((x2+40.021433)*x2+322.62491)*x2+570.23628)*x2+ ...
              157.105423);
        gx = ((((x2+42.242855)*x2+302.757865)*x2+352.018498)*x2+21.821899)/ ...
             ((((x2+48.196927)*x2+482.485984)*x2+1114.978885)*x2+ ...
              449.690326)/x;
        y = 1.570796327-fx*cos(x)/x-gx*sin(x)/x;
    end
end
