function y = fftconv(x,h,mode)
%y = fftconv(x,h,mode)
%fft convolution is much faster than overlap-add
%   x,h the signals to be convolved
%   mode = 1 zero pads h at the tail end (default)
%   mode = 2 zero pads h symmetrically

if exist('mode') ~= 1
    mode = 1;
end

lx = length(x);
lh = length(h);

if mode == 2
    y = real(ifft(fft([x zeros(1,lh-1)]).*fft([zeros(1,floor(lx/2)) h zeros(1,ceil(lx/2)-1)])));
else
    y = real(ifft(fft([x zeros(1,lh-1)]).*fft([h zeros(1,lx-1)])));
end

%if any(isnan(y))
%    'Reverting to manual convolution... this might take a while.'
%    y = conv(x,h);
%end