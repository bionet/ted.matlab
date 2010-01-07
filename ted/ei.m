%EI Exponential integral function.
%   Y = EI(Z) computes the exponential integral Ei of 
%   a complex value.

%   Adapted from specfun (http://www.netlib.org/specfun/ 
%   by Lev Givon
%   Copyright 2009 Trustees of Columbia University


function y = ei(z)

if any(not(isreal(z))),
    y = -exp1(-z)+(log(z)-log(1./z))/2-log(-z);
else
    y = real(-exp1(-z)+(log(z)-log(1./z))/2-log(-z));
end
