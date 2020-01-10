clear; close all
load('tmpdata')
x = dat(1,:);

dx = fft(x);

N = size(dx,2);
dx_1s = dx(1:(N/2)+1);
dx_2s = [dx_1s fliplr(conj(dx_1s(:,2:end-1)))];

plot(abs(dx_2s)); hold on; plot(abs(dx))
figure
plot(angle(dx_2s)); hold on; plot(angle(dx))
