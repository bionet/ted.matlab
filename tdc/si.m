%SI Sine integral function.
%   Y = SI(X) computes the sin integral INT(SIN(T)/T,T,0,X)
%   using Bessel functions. This implementation is faster than
%   Matlab's SININT function.

%   Adapted from Si.m by Sylvain Pelissier <sylvain.pelissier@gmail.com>
%   This file is licensed under GPL version 2 or later.

function y = si(x)

y = zeros(size(x));
n = prod(size(x));
if n < 101,
  for k=1:n,
    y(k) = sum(besselj([0:100]+0.5,(x(k)/2)).^2);
  end
else
  for k=0:100,
    y = y + besselj(k+0.5,x/2).^2;
  end
end 
y = real(y).*pi;

