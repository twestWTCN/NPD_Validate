function freqS = freqShuff(freq)
freqS = freq;
freqshuff = freq.fourierspctrm;
for i = 1:size(freqshuff,2)
    X = squeeze(freqshuff(:,i,:));
    freqshuff(:,i,:) = X(randperm(size(X,1)),:);
end
freqS.fourierspctrm = freqshuff;
