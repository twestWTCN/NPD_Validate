function freqS = phaseRandSym(freq)

freqshuff = freq;
for i = 1:size(freqshuff,2)
    X = squeeze(freqshuff(i,:));
    if rem(size(X,2),2)
        Y = X(1,1:floor(size(X,2)/2));
    
    for j = 1:size(X,1)
        Y = X(j,:);
        rp = 2*pi*rand(1,size(X,2));
        X(j,:) = Y.*exp(sqrt(-1).*rp);
    end
    freqshuff(:,i,:) = X;
end
freqS.fourierspctrm = freqshuff;
