function freqS = freqShuff(freq,bstraptype)
freqS = freq;
freqshuff = freq.fourierspctrm;

if bstraptype == 1
    % shuffle all
    for i = 1:size(freqshuff,2)
        iShuff = randi(size(freqshuff,1),1,size(freqshuff,1));
        freqshuff(:,i,:) = freqshuff(iShuff,i,:);
    end
end

freqS.fourierspctrm = freqshuff;
