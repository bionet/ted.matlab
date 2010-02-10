function [h,fc,t,f] = gammatone(num,len,fmin,fmax,fs,pad_bw)
%[h,fc,t,f] = gammatone(num,len,fmin,fmax,fs,pad_bw)
%Create gammatone filterbank
%   num = number of filters
%   len = length, in samples, of the filter
%   fmin = center frequency of first filter
%   fmax = maximum frequency of the filterbank
%          fc(end) = fmax if pad_bw = 1
%   fs   = sampling frequency
%   pad_bw = Padding of the highest band to avoid aliasing. Default=2.
%Returns:
%   h  = numxlen array of filters
%   fc = vector of center frequencies
%   t  = time used to plot filterbank
%   f  = frequency used to plot filterbank

if exist('pad_bw') ~= 1
    pad_bw = 2;
end

EarQ = 9.26449;
minBW = 24.7;
order = 4;
dt = 1/fs;
t = dt*(0:len-1);
f = (0:length(t)-1)/length(t)*fs;
beta = 1.019;

Wp = fmax;
fmax = EarQ*(Wp-pad_bw*beta*minBW)/(EarQ+pad_bw*beta);

overlap = EarQ*(log(fmax+EarQ*minBW)-log(fmin+EarQ*minBW))/max(1,num-1);
fc = -EarQ*minBW + (fmax+EarQ*minBW)*exp(-(num-(1:num))*overlap/EarQ)
h = zeros(num,len);
for i=1:num
    h(i,:) = t.^(order-1).*exp(-2*pi*beta*(fc(i)/EarQ+minBW)*t).*cos(2*pi*fc(i)*t);
    h(i,:) = h(i,:)/max(abs(fft(h(i,:))));  %alpha
end
