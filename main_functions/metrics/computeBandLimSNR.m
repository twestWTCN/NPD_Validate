function snr = computeBandLimSNR(s,n,bdef,dataform)
dataf.label = {'sig','noise'};
dataf.trial{1} = vertcat(s,n);
dataf.time = dataform.time;
dataf.fsample = dataform.fsample;
spec = computeSpectra(dataf,[],2,0,[],-1,8);%(dataN,[0 0 0],Nsig,plotfig,linestyle,-1,datalength)
inds = find(spec.freq>=bdef(1) & spec.freq<=bdef(2));

snr = sum(mean(abs(squeeze(spec.fourierspctrm(:,1,inds))).^2))./sum(mean(abs(squeeze(spec.fourierspctrm(:,2,inds)).^2)));
