X = data_macro.trial{1}(1,:);
fsamp = data_macro.fsample;
periodogram(X,[],fsamp,fsamp);
hold on
Xshuff  = X(randperm(size(X,2)));
periodogram(Xshuff,[],fsamp,fsamp);

Xw = fft(X);
Xwshuff = ifft(Xw(randperm(size(Xw,2))));
periodogram(Xwshuff,[],fsamp,fsamp);

figure
subplot(2,1,1)
plot(X)
hold on
plot(Xshuff)
plot(Xwshuff)

std(X)
std(Xshuff)
std(Xwshuff)