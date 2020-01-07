function freqS = phaseRand(freq)
freqS = freq;
freqshuff = freq.fourierspctrm;
for i = 1:size(freqshuff,2)
    X = squeeze(freqshuff(:,i,:));
    for j = 1:size(X,1)
        Y = X(j,:);
        rp = 2*pi*rand(1,size(X,2));
        X(j,:) = Y.*exp(sqrt(-1).*rp);
    end
    freqshuff(:,i,:) = X;
end
freqS.fourierspctrm = freqshuff;
